/* This header file contains functions common to both qfind and qfind-s.
** Some functions contained in this file work slightly differently depending
** on whether qfind or qfind-s is being compiled.  Such differences are
** determined by the presence of the macro QSIMPLE defined in qfind-s.c.
*/

#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE   /* MAP_ANONYMOUS / madvise under -std=c11 */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <stdatomic.h>
#ifdef __linux__
#include <sys/mman.h>
#endif

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
   #define omp_get_num_procs() 2
   #define omp_set_num_threads(x)
#endif

#ifdef QSIMPLE
   #define WHICHPROGRAM qfind-simple
#else
   #define WHICHPROGRAM qfind
#endif

#define STR(x) #x
#define XSTR(x) STR(x)
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define BANNER XSTR(WHICHPROGRAM)" v2.4b by Matthias Merzenich, 8 September 2025"

#define FILEVERSION ((unsigned long) 2025090401)  /* yyyymmddnn */

#define MAXPERIOD 30
#define MAXDUMPROOT 50     /* maximum allowed length of dump root */
#define DUMPLIMIT 100000   /* maximum allowed number of sequential dumps */
#define CHUNK_SIZE 1
#define QBITS 15
#define HASHBITS 15
#define DEFAULT_DEPTHLIMIT (qBits-3)
#define DEFAULT_CACHEMEM 32

#define P_WIDTH 0
#define P_PERIOD 1
#define P_OFFSET 2
#define P_SYMMETRY 3
#define P_REORDER 4           /* currently cannot be set at runtime */
#define P_DUMPMODE 5
#define P_BASEBITS 6
#define P_QBITS 7
#define P_HASHBITS 8
#define P_DEPTHLIMIT 9        /* currently cannot be set at runtime */
#define P_NUMTHREADS 10
#define P_MINDEEP 11
#define P_MEMLIMIT 12
#define P_CACHEMEM 13
#define P_PRINTDEEP 14
#define P_LONGEST 15
#define P_FIRSTDEEP 16
#define P_NUMSHIPS 17
#define P_MINEXTENSION 18
#define P_FULLPERIOD 19
#define P_BOUNDARYSYM 20
#define P_DUMPINTERVAL 21
#define P_EVERYDEPTH 22
#define P_EARLYEXIT 23

#define NUM_PARAMS 24U

#define SYM_UNDEF 0
#define SYM_ASYM 1
#define SYM_ODD 2
#define SYM_EVEN 3
#define SYM_GUTTER 4

const char *rule = "B3/S23";     /* Default rule set to B3/S23 (Life) */
char loadRule[151];              /* Used for loading rule from file */
char baseRule[151];              /* Used in case of forbidden conditions */

char *initRows;

int params[NUM_PARAMS];
int width;
int nRowsInState;                /* Could be replaced with 2*period.  Should be removed. */

int period;
int offset;
int lastDeep = 0;
int numFound = 0;       /* number of spaceships found so far */
int longest = 0;        /* length of current longest partial result */

int aborting = 0;       /* A flag to indicate that we are stopping the search.
                        ** valid aborting values:
                        **   0: not aborting
                        **   1: fatal error
                        **   2: queue size limit reached
                        **   3: desired number of ships found
                        */

int gutterSkew = 0;
int lookaheadDepth = 3;  /* --lookahead-depth: 3 = full lookahead; 1-2 = weaker, cheaper */     /* number of cells to skew halves in gutter symmetric search */

/* the big data structures */
#define qBits params[P_QBITS]
#define QSIZE (1LLU<<qBits)

#define hashBits params[P_HASHBITS]
#define HASHSIZE (1LU<<hashBits)
#define HASHMASK (HASHSIZE - 1)

typedef uint32_t node;   /* index of a node (partial pattern) in the BFS queue */
typedef uint16_t row;    /* one row of cells, one bit per cell (width <= 14)   */

row * queueRows;         /* per queue node: its last row + offset to its parent (see below) */
node * queueBase;        /* per block of 2^BASEBITS nodes: base index the parent offsets add to */
node * seenHash;         /* hash of recently seen states (last 2*period rows) for duplicate removal */

int ruleTable[512];      /* rule: 9-cell neighbourhood -> next state (1, 0, or -1 = forbidden) */
_Atomic(uint16_t*) *pairTable;   /* lookup table per pair of rows (row1,row2), built on first use */
uint32_t *rowOrderWeight;        /* heuristic weight per row; candidate rows are tried heaviest first */
long long memusage = 0;          /* memory accounted by qfind itself (tables, index, hash, cache) */
long long memlimit = 0;          /* -m limit in bytes (only used when -m >= 0) */

/* Saved depth-first extensions.  When the deepening step finds that a queue node can
** be extended, the rows it found are saved (extensionRows[i] = {last index, start
** index, rows...}) so the next deepening or BFS step can resume from them instead of
** searching again.  extensionIndex runs parallel to the queue (from extQHead): 0 = no
** extension, 1 = extension found but not saved, >1 = index into extensionRows.     */
row **extensionRows = 0;
/* extensionIndex needs 2*QSIZE + 1 entries: while doCompactPart2() re-enqueues the   */
/* whole queue, enqueue() keeps advancing extQTail from its old value (up to QSIZE)  */
/* and clears the entry there.  (The original allocated 2*QSIZE entries by accident, */
/* via sizeof(pointer), which hid this.)                                             */
#define EXTINDEX_SIZE (2*QSIZE + 1)
uint32_t *extensionIndex;
uint32_t extQHead, extQTail, prevExtQHead;

#ifndef NOCACHE
long long cachesize;
typedef struct {
   uint64_t id1, id2, id3;   /* identities of list11, list12, list13 */
   int rowsKey;              /* the two remaining rows the result depends on */
   int r;                    /* cached lookahead result (0 or 1) */
} cacheentry;

cacheentry *totalCache;

cacheentry **cache;
#endif

/* Representation of vertices.
**
** Each vertex (node) is represented by an entry in the queueRows[] array.
** That entry consists of the bits of the last row of the vertex's pattern
** concatenated (above bit `width`) with a number, the "offset" of the vertex.
** The parent of vertex i is formed by adding the offset of i
** to the number in queueBase[i/BASEFACTOR].  A node's pattern is recovered by
** following parents back to the root (node 0).
**
** If node i has offset -1 (all bits set), the entry is considered to be empty and
** is skipped.  This is to deal with the case that queueBase[i/BASEFACTOR] is already
** set to something too far away to deal with via the offset.
**
** qIsEmpty() checks the size of the queue
** enqueue(n,r) adds a new node with last row r and parent n
** dequeue() returns the index of the next unprocessed node
** pop() undoes the previous enqueue operation
** resetQ() clears the queue
**
** ROW(b) returns the last row of b
** PARENT(b) returns the index of the parent of b
*/

#define MAXWIDTH (14)   /* hard limit as long as rows are of type uint16_t */

#define ROWBITS ((1<<width)-1)
#define BASEBITS (params[P_BASEBITS])
#define BASEFACTOR (1<<BASEBITS)
#define MAX_PARENT_OFFSET ((((row) -1) >> width) - 1)

#define ROW(i) (queueRows[i] & ROWBITS)
#define PARENT_OFFSET(i) (queueRows[i] >> width)
#define EMPTY(i) (queueRows[i] == (row)-1)
#define MAKEEMPTY(i) queueRows[i] = (row)-1
#define PARENT(i) (queueBase[(i)>>BASEBITS]+PARENT_OFFSET(i))
#define FIRSTBASE(i) (((i) & ((1<<BASEBITS) - 1)) == 0)

#define MINDEEP ((params[P_MINDEEP]>0) ? params[P_MINDEEP] : 3)

void optError(const char *errorMsg, const char *opt);
void printError(const char *errorMsg);

int gcd(int a, int b) {
   if (a > b) return gcd(b,a);
   else if (a == 0) return b;
   else return gcd(b-a,a);
}

/* =============================== */
/*  Display current date and time  */
/* =============================== */

static char timeStr[20] = "00/00/00 00:00:00 ";

static void timeStamp(void) {
   time_t t;
   
   time(&t);
   strftime(timeStr, 20, "%d/%m/%y %H:%M:%S ", localtime(&t));
   printf("%s", timeStr);
}


/* ---- progress readout: elapsed wall time and memory ---- */
#include <sys/resource.h>
#include <unistd.h>
static struct timespec startClock;
static void __attribute__((constructor)) startClockInit(void) { clock_gettime(CLOCK_MONOTONIC, &startClock); }
static double elapsedSeconds(void) {
   struct timespec n;
   clock_gettime(CLOCK_MONOTONIC, &n);
   return (n.tv_sec - startClock.tv_sec) + (n.tv_nsec - startClock.tv_nsec) / 1e9;
}
/* resident set size now (Linux) in MB, or -1 if unknown */
static double currentRssMB(void) {
#ifdef __linux__
   FILE *f = fopen("/proc/self/statm", "r");
   long pages = 0, res = 0;
   if (f) { if (fscanf(f, "%ld %ld", &pages, &res) != 2) res = 0; fclose(f); }
   return res ? res * (double)sysconf(_SC_PAGESIZE) / 1048576.0 : -1;
#else
   return -1;
#endif
}
/* peak resident set size in MB */
static double peakRssMB(void) {
   struct rusage ru;
   getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
   return ru.ru_maxrss / 1048576.0;   /* bytes */
#else
   return ru.ru_maxrss / 1024.0;      /* kilobytes */
#endif
}
static void printProgress(void) {
   double cur = currentRssMB();
   if (cur >= 0) printf("  [%.1fs, mem %.0f MB (peak %.0f), tables+index %lld MB]", elapsedSeconds(), cur, peakRssMB(), memusage >> 20);
   else          printf("  [%.1fs, peak mem %.0f MB, tables+index %lld MB]", elapsedSeconds(), peakRssMB(), memusage >> 20);
}

/* =========================================== */
/*  Lookup tables to determine successor rows  */
/* =========================================== */

const char *rulekeys[] = {
   "0", "1c", "1e", "2a", "1c", "2c", "2a", "3i",
   "1e", "2k", "2e", "3j", "2a", "3n", "3a", "4a",
   "1c", "2n", "2k", "3q", "2c", "3c", "3n", "4n",
   "2a", "3q", "3j", "4w", "3i", "4n", "4a", "5a",
   "1e", "2k", "2i", "3r", "2k", "3y", "3r", "4t",
   "2e", "3k", "3e", "4j", "3j", "4k", "4r", "5n",
   "2a", "3q", "3r", "4z", "3n", "4y", "4i", "5r",
   "3a", "4q", "4r", "5q", "4a", "5j", "5i", "6a",
   "1c", "2c", "2k", "3n", "2n", "3c", "3q", "4n",
   "2k", "3y", "3k", "4k", "3q", "4y", "4q", "5j",
   "2c", "3c", "3y", "4y", "3c", "4c", "4y", "5e",
   "3n", "4y", "4k", "5k", "4n", "5e", "5j", "6e",
   "2a", "3n", "3r", "4i", "3q", "4y", "4z", "5r",
   "3j", "4k", "4j", "5y", "4w", "5k", "5q", "6k",
   "3i", "4n", "4t", "5r", "4n", "5e", "5r", "6i",
   "4a", "5j", "5n", "6k", "5a", "6e", "6a", "7e",
   "1e", "2a", "2e", "3a", "2k", "3n", "3j", "4a",
   "2i", "3r", "3e", "4r", "3r", "4i", "4r", "5i",
   "2k", "3q", "3k", "4q", "3y", "4y", "4k", "5j",
   "3r", "4z", "4j", "5q", "4t", "5r", "5n", "6a",
   "2e", "3j", "3e", "4r", "3k", "4k", "4j", "5n",
   "3e", "4j", "4e", "5c", "4j", "5y", "5c", "6c",
   "3j", "4w", "4j", "5q", "4k", "5k", "5y", "6k",
   "4r", "5q", "5c", "6n", "5n", "6k", "6c", "7c",
   "2a", "3i", "3j", "4a", "3q", "4n", "4w", "5a",
   "3r", "4t", "4j", "5n", "4z", "5r", "5q", "6a",
   "3n", "4n", "4k", "5j", "4y", "5e", "5k", "6e",
   "4i", "5r", "5y", "6k", "5r", "6i", "6k", "7e",
   "3a", "4a", "4r", "5i", "4q", "5j", "5q", "6a",
   "4r", "5n", "5c", "6c", "5q", "6k", "6n", "7c",
   "4a", "5a", "5n", "6a", "5j", "6e", "6k", "7e",
   "5i", "6a", "6c", "7c", "6a", "7e", "7c", "8" };

/*   Parses the rule.  If there's an error, return a string describing the error.
**   Fills in the 512-element array pointed to by tab.
*/
const char *parseRule(const char *theRule, int *table) {
   const char *p = theRule;
   int letterHits[256];    /* needed to keep track of values when setting forbidden conditions */
   for (int i=0; i<512; i++)
      table[i] = 0;
   for (int ruleOffset=0; ruleOffset<512; ruleOffset += 256) {
      if (ruleOffset == 0) {
         if (*p != 'B' && *p != 'b')
            return "Expected B at start of rule";
      } else {
         if (*p != 'S' && *p != 's')
            return "Expected S after slash";
      }
      p++;
      int allowed = 1;
      while (*p != '/' && *p != '\0') {
         if (*p == '~'){
            p++;
            if (allowed == -1 || *p == '~'){
               if (ruleOffset)
                  return "Can't have multiple tildes in survival conditions";
               else
                  return "Can't have multiple tildes in birth conditions";
            }
            if (*p == '/' || *p == '\0')
               continue;
            allowed = -1;  /* table entry will be -1 if condition is forbidden */
         }
         if (!('0' <= *p && *p <= '9'))
            return "Missing number in rule";
         if (*p == '9')
            return "Unexpected character in rule";
         char digit = *p++;
         int negated = 0;
         if (*p == '/' || *p == '\0' || (*p == '-' && allowed == 1) || *p == '~' || ('0' <= *p && *p <= '8'))
            for (int i=0; i<256; i++)
               if (rulekeys[i][0] == digit)
                  table[ruleOffset+i] = 1*allowed;
         int forbiddenCount = 0;
         if (*p == '-'){
            negated = 1;
            for (int i=0; i<256; i++)
               letterHits[i] = 0;
            p++;
         }
         for (; *p != '/' && *p != '\0' && *p != '~' && !('0' <= *p && *p <= '8'); p++) {
            if (*p == '-')
               return "Improperly placed negation sign";
            if ('a' <= *p && *p <= 'z') {
               int used = 0;
               for (int i=0; i<256; i++){
                  if (rulekeys[i][0] == digit){
                     if (rulekeys[i][1] == *p){
                        if (allowed == 1)
                           table[ruleOffset+i] = (1-negated);
                        else if (!negated)
                           table[ruleOffset+i] = -1;
                        used++;
                     }
                     else if (negated && allowed == -1)
                        letterHits[i]++;
                  }
               }
               if (negated && allowed == -1)
                  forbiddenCount++;
               if (!used)
                  return "Unexpected character in rule";
            }
            else
               return "Unexpected character in rule";
         }
         if (negated && (*(p-1) == '-'))
            return "Negation sign must be followed by at least one letter";
         if (negated && allowed == -1)
            for (int i=0; i<256; i++)
               if (letterHits[i] == forbiddenCount)
                  table[ruleOffset+i] = -1;
      }
      if (ruleOffset == 0) {
         if (*p++ != '/')
            return "Missing expected slash between B and S";
      } else {
         if (*p++ != 0)
            return "Extra unparsed junk at end of rule string";
      }
   }
   return 0;
}

#ifndef QSIMPLE
void makePhases(void);
#endif

int8_t ruleTable3x3[512];

int slowEvolveBit(int row1, int row2, int row3, int bshift) {
   return ruleTable[(((row2>>bshift) & 2)<<7) | (((row1>>bshift) & 2)<<6)
                | (((row1>>bshift) & 4)<<4) | (((row2>>bshift) & 4)<<3)
                | (((row3>>bshift) & 7)<<2) | (((row2>>bshift) & 1)<<1)
                |  ((row1>>bshift) & 1)<<0];
}

void makeRuleTable3x3(void) {
   int p = 0;
   for (int row1=0; row1<8; row1++)
      for (int row2=0; row2<8; row2++)
         for (int row3=0; row3<8; row3++)
            ruleTable3x3[p++] = (int8_t)slowEvolveBit(row1, row2, row3, 0);
}

int evolveBitShift(int row1, int row2, int row3, int bshift) {
   return ruleTable3x3[
      (((row1 << 6) >> bshift) & 0700) +
      (((row2 << 3) >> bshift) &  070) +
      (( row3       >> bshift) &   07)];
}

int evolveBit(int row1, int row2, int row3) {
   return ruleTable3x3[
      ((row1 << 6) & 0700) +
      ((row2 << 3) &  070) +
      ( row3       &   07)];
}

int evolveRow(int row1, int row2, int row3) {
   int row4;
   int row1Shifted,row2Shifted,row3Shifted;
   int j,symShift = 0;
   int boundaryShift = 0;
   int theBit;
   if (params[P_BOUNDARYSYM] == SYM_GUTTER && !gutterSkew){ /* Need to check the gutter for forbidden births */
      theBit = (row1 >> (width-1)) + ((row2 >> (width-1)) << 1) + ((row3 >> (width-1)) << 2);
      if (evolveBit(theBit, 0, theBit))
         return -1;
   }
   if (params[P_SYMMETRY] == SYM_GUTTER && !gutterSkew){    /* Need to check the gutter for forbidden births */
      theBit = (row1 & 1) + ((row2 & 1) << 1) + ((row3 & 1) << 2);
      if (evolveBit(theBit, 0, theBit))
         return -1;
   }
   if (params[P_SYMMETRY] == SYM_ODD) symShift = 1;
   if (params[P_BOUNDARYSYM] == SYM_UNDEF && evolveBitShift(row1, row2, row3, width - 1)) return -1;
   if (params[P_BOUNDARYSYM] == SYM_ODD) boundaryShift = 1;
   if (params[P_SYMMETRY] == SYM_ASYM && evolveBit(row1 << 2, row2 << 2, row3 << 2)) return -1;
   if (params[P_SYMMETRY] == SYM_ODD || params[P_SYMMETRY] == SYM_EVEN){
      row1Shifted = (row1 << 1) + ((row1 >> symShift) & 1);
      row2Shifted = (row2 << 1) + ((row2 >> symShift) & 1);
      row3Shifted = (row3 << 1) + ((row3 >> symShift) & 1);
   }
   else {
      row1Shifted = (row1 << 1);
      row2Shifted = (row2 << 1);
      row3Shifted = (row3 << 1);
   }
   if (params[P_BOUNDARYSYM] == SYM_ODD || params[P_BOUNDARYSYM] == SYM_EVEN){
      row1 += ((row1 >> (width-1-boundaryShift)) & 1) << (width);
      row2 += ((row2 >> (width-1-boundaryShift)) & 1) << (width);
      row3 += ((row3 >> (width-1-boundaryShift)) & 1) << (width);
   }
   row4 = evolveBit(row1Shifted, row2Shifted, row3Shifted);
   if (row4 == -1) return -1;
   for (j = 1; j < width; j++){
      theBit = evolveBitShift(row1, row2, row3, j - 1);
      if (theBit == -1) return -1;
      row4 += theBit << j;
   }
   return row4;
}

int evolveRowHigh(int row1, int row2, int row3, int bits) {
   int row4=0;
   int j,boundaryShift = 0;
   int theBit;
   if (params[P_BOUNDARYSYM] == SYM_GUTTER && !gutterSkew){ /* Need to check the gutter for forbidden births */
      theBit = (row1 >> (width-1)) + ((row2 >> (width-1)) << 1) + ((row3 >> (width-1)) << 2);
      if (evolveBit(theBit, 0, theBit))
         return -1;
   }
   if (params[P_BOUNDARYSYM] == SYM_UNDEF && evolveBitShift(row1, row2, row3, width - 1)) return -1;
   if (params[P_BOUNDARYSYM] == SYM_ODD) boundaryShift = 1;
   if (params[P_BOUNDARYSYM] == SYM_ODD || params[P_BOUNDARYSYM] == SYM_EVEN){
      row1 += ((row1 >> (width-1-boundaryShift)) & 1) << (width);
      row2 += ((row2 >> (width-1-boundaryShift)) & 1) << (width);
      row3 += ((row3 >> (width-1-boundaryShift)) & 1) << (width);
   }
   for (j = width-bits; j < width; j++){
      theBit = evolveBitShift(row1, row2, row3, j - 1);
      if (theBit == -1)
         return -1;
      row4 += theBit << j;
   }
   return row4;
}

int evolveRowLow(int row1, int row2, int row3, int bits) {
   int row4;
   int row1Shifted,row2Shifted,row3Shifted;
   int j,symShift = 0;
   int theBit;
   if (params[P_SYMMETRY] == SYM_GUTTER && !gutterSkew){ /* Need to check the gutter for forbidden births */
      theBit = (row1 & 1) + ((row2 & 1) << 1) + ((row3 & 1) << 2);
      if (evolveBit(theBit, 0, theBit))
         return -1;
   }
   if (params[P_SYMMETRY] == SYM_ODD) symShift = 1;
   if (params[P_SYMMETRY] == SYM_ASYM && evolveBit(row1 << 2, row2 << 2, row3 << 2)) return -1;
   if (params[P_SYMMETRY] == SYM_ODD || params[P_SYMMETRY] == SYM_EVEN){
      row1Shifted = (row1 << 1) + ((row1 >> symShift) & 1);
      row2Shifted = (row2 << 1) + ((row2 >> symShift) & 1);
      row3Shifted = (row3 << 1) + ((row3 >> symShift) & 1);
   }
   else{
      row1Shifted = (row1 << 1);
      row2Shifted = (row2 << 1);
      row3Shifted = (row3 << 1);
   }
   row4 = evolveBit(row1Shifted, row2Shifted, row3Shifted);
   if (row4 == -1) return -1;
   for (j = 1; j < bits; j++){
      theBit = evolveBitShift(row1, row2, row3, j - 1);
      if (theBit == -1) return -1;
      row4 += theBit << j;
   }
   return row4;
}

/* sortRows() is used to sort the global array rowOrder. 
** This array determines the order in which new rows are 
** added in the search.  The idea is that certain search
** orders give slightly quicker solutions during depthFirst()
** and LookAhead() than the naive order.
*/
void sortRows(uint16_t *theRow, uint32_t totalRows) {
   uint32_t i;
   int64_t j;
   uint16_t t;
   for (i = 1; i < totalRows; ++i){
      t = theRow[i];
      j = i - 1;
      while (j >= 0 && rowOrderWeight[theRow[j]] < rowOrderWeight[t]){
         theRow[j+1] = theRow[j];
         --j;
      }
      theRow[j+1] = t;
   }
}

/* Lookup table layout (16-bit words):
**
**    mask[MASKWORDS]   bit r set iff bucket r (rows row3 that produce row r) is non-empty
**    pre[PREWORDS]     number of non-empty buckets before each 64-bit chunk of the mask
**    total             total number of words in the table
**    dbase             word index where the list data begins
**    off[ne+1]         start of every NON-EMPTY bucket's list (and the end of the last one)
**    data              the row lists, in bucket order
**
** Offsets are stored only for the non-empty buckets and found with a rank (popcount)
** query.  With PACKED (opt-in: compile with -DPACKLISTS) every list
** entry is stored in `width` bits instead of 16 and `off` counts entries instead of
** words.  A list is then decoded into a caller-supplied buffer of DECODE_CAP entries
** (bucketEntries).  A table pointer points at mask[0].                                       */
#ifdef PACKLISTS
   #define PACKED 1
#endif
#define MASKWORDS ((((1<<width) + 63) >> 6) << 2)
#define PREWORDS (MASKWORDS >> 2)
#define HDRWORDS (MASKWORDS + PREWORDS)
#define DECODE_CAP ((1<<width) + 8)
/* words of per-thread scratch for building one table: header (HDRWORDS + 2), at most  */
/* 2^width + 1 offsets, at most 2^width list entries, plus 16 words of padding.     */
#define SCRATCH_WORDS ((2LL << width) + HDRWORDS + 19)
#define BITSET(m, r) (((m)[(r) >> 4] >> ((r) & 15)) & 1)

/* index (among the non-empty buckets) of non-empty bucket r */
static inline int bucketRank(const uint16_t *t, int r) {
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
   uint64_t chunk;
   memcpy(&chunk, t + ((r >> 6) << 2), sizeof chunk);
   return t[MASKWORDS + (r >> 6)] + __builtin_popcountll(chunk & ((1ULL << (r & 63)) - 1));
#else
   int idx = t[MASKWORDS + (r >> 6)];
   for (int w = (r >> 6) << 2; w < (r >> 4); w++) idx += __builtin_popcount(t[w]);
   return idx + __builtin_popcount(t[r >> 4] & ((1u << (r & 15)) - 1));
#endif
}

/* bucket r of table t: returns its length, *start = position of its first entry */
static inline int bucketFind(const uint16_t *t, int r, int *start) {
   if (!BITSET(t, r)) { *start = 0; return 0; }
   const uint16_t *off = t + HDRWORDS + 2;
   int idx = bucketRank(t, r);
   *start = off[idx];
   return off[idx + 1] - off[idx];
}

/* the entries start..start+n-1 of table t, written to / pointed to via buf */
static inline const uint16_t *bucketRead(const uint16_t *t, int start, int n, uint16_t *buf) {
#ifdef PACKED
   const uint16_t *d = t + t[HDRWORDS + 1];
   const int w = width;
   const uint32_t m = (1u << w) - 1;
   unsigned bp = (unsigned)start * w;
   int i = 0;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
   /* 4 entries per 64-bit load (7 bits of misalignment + 4*14 bits <= 64).  Always a
   ** whole number of groups is written, i.e. up to 3 entries past the end of the list
   ** (buffers have DECODE_CAP >= 2^width + 8 entries; the table data is padded).        */
   {
      uint64_t v;
      memcpy(&v, (const uint8_t*)d + (bp >> 3), 8);
      v >>= (bp & 7);
      buf[0] = (uint16_t)(v & m);
      buf[1] = (uint16_t)((v >> w) & m);
      buf[2] = (uint16_t)((v >> (2 * w)) & m);
      buf[3] = (uint16_t)((v >> (3 * w)) & m);
      i = 4; bp += 4 * w;
   }
   for (int g = (n - 1) >> 2; g > 0; g--, i += 4, bp += 4 * w) {
      uint64_t v;
      memcpy(&v, (const uint8_t*)d + (bp >> 3), 8);
      v >>= (bp & 7);
      buf[i]     = (uint16_t)(v & m);
      buf[i + 1] = (uint16_t)((v >> w) & m);
      buf[i + 2] = (uint16_t)((v >> (2 * w)) & m);
      buf[i + 3] = (uint16_t)((v >> (3 * w)) & m);
   }
#else
   for (; i < n; i++, bp += w) {
      unsigned wi = bp >> 4;
      uint32_t v = d[wi] | ((uint32_t)d[wi + 1] << 16);
      buf[i] = (uint16_t)((v >> (bp & 15)) & m);
   }
#endif
   return buf;
#else
   (void)buf; (void)n;
   return t + start;
#endif
}

/* the list of row3 values in bucket r of table t: sets *n; buf needs DECODE_CAP entries */
static inline const uint16_t *bucketEntries(const uint16_t *t, int r, int *n, uint16_t *buf) {
   int start;
   *n = bucketFind(t, r, &start);
   if (!*n) return buf;
   return bucketRead(t, start, *n, buf);
}

/* total size in words of a stored table */
static inline int tableWords(const uint16_t *t) { return t[HDRWORDS]; }

uint16_t *buildTable(int row1, int row2);

/* tableForPair() returns a pointer to a lookup table where further information */
/* is found.  It is used in getSuccessors() and in lookAhead().             */
uint16_t *tableForPair(int row12) {
   /* acquire pairs with the release store in buildTable(): a thread that sees the  */
   /* pointer also sees the table contents (matters on ARM, e.g. Apple Silicon).    */
   uint16_t *r = atomic_load_explicit(&pairTable[row12], memory_order_acquire);
   if (r == 0)
      r = buildTable(row12 >> width, row12 & ((1 << width) - 1));
   return r;
}
uint16_t *tableFor(int row1, int row2) {
   return tableForPair((row1 << width) + row2);
}

/* ========================================================================== *
** Two representations of the lookup tables; each search uses one of them,
** chosen at startup (--tables, default: per-pair below width 10, factored from 10).
**
**  per-pair (the original idea): one table per pair (row1,row2), built on first
**      use (buildTable) and kept.  Fastest while the tables fit in fast memory, but
**      a real search builds nearly all 4^width of them: about 1.5 GB at width 10,
**      roughly 7x more per extra width.
**
**  factored: the rule is evaluated in two halves that overlap by 2 bits.  The low
**      result bits (0 .. lowResultBits-1) of row4 depend only on the low lowHalfBits bits of row1,
**      row2 and row3; the high result bits (lowResultBits .. width-1) only on bits
**      highHalfShift .. width-1 of the three rows (highHalfShift = lowHalfBits - 2, so row3 bits highHalfShift and highHalfShift+1
**      belong to both halves).  buildTable() already used exactly this split; here
**      the two halves are stored globally, once, instead of combining them into a
**      table per row pair:
**        low:  for each (row1 low, row2 low): the row3-low values grouped by
**              (row4 low bits, overlap value v)     -> lowHalfOffsets / lowHalfRows, lowHalfMasks
**        high: the same for the high halves          -> highHalfOffsets / highHalfRows, highHalfMasks
**      A successor list for (row1, row2, row4) is the join of the two groups with
**      equal overlap value; "is the list non-empty" is one AND of two 4-bit masks.
**      Memory grows about 2.8x per width instead of about 7x: 3 MB at width 10,
**      about 170 MB at width 14.  Lists come out in a different order than the
**      per-pair tables (rows of one list are the same); that changes timing and which
**      variants of a ship get printed, never whether a ship is found.  Compile with
**      -DFACTORED_SORTED to sort every list into the per-pair order (slower; used to
**      check that both representations give identical searches).
**
** A TableRef is a table pointer (per-pair) or a packed pair of half-table indices
** (factored).  TABLE_REF / TABLE_HAS_ROW / TABLE_ROWS take the representation as their first
** argument; the hot lookahead is compiled once for each so neither pays a branch.
** ========================================================================== */
typedef uintptr_t TableRef;
int useFactored = 0;         /* representation chosen for this search */
int tableMode = 0;           /* --tables: 0 = auto, 1 = per-pair, 2 = factored */
#ifndef LOOKAHEAD_BUDGET
#define LOOKAHEAD_BUDGET 256
#endif
int lookaheadMethod = 0;     /* --lookahead-method: 0 = auto, 1 = cone, 2 = split, 3 = hybrid */
int useSplitLookahead = 0;   /* chosen for this search: 0 = cone, 1 = split, 2 = hybrid (1, 2 need the half-tables) */
long lookaheadBudget = LOOKAHEAD_BUDGET;   /* --lookahead-budget: work units before the hybrid switches to split */
#ifndef SPLIT_CACHE_BITS
#define SPLIT_CACHE_BITS 22
#endif
int splitCacheBits = SPLIT_CACHE_BITS;     /* --split-cache-bits: 2^N shared split-cache entries per half */
static int lowHalfBits, lowResultBits, highHalfShift, highHalfBits, highResultBits;
static uint8_t *lowHalfRows, *highHalfRows, *lowHalfMasks, *highHalfMasks;
static uint16_t *lowHalfOffsets, *highHalfOffsets;
static uint32_t *pairTableRank;
extern uint16_t *rowOrder;
#define LOW_OFFSETS_PER_PAIR ((4 << lowResultBits) + 1)
#define HIGH_OFFSETS_PER_PAIR ((4 << highResultBits) + 1)

static inline uint32_t factoredTableRef(int row1, int row2) {
   const int lowMask = (1 << lowHalfBits) - 1;
   uint32_t lowPair = ((uint32_t)(row1 & lowMask) << lowHalfBits) | (uint32_t)(row2 & lowMask);
   uint32_t highPair = ((uint32_t)(row1 >> highHalfShift) << highHalfBits) | (uint32_t)(row2 >> highHalfShift);
   return lowPair | (highPair << 16);
}
static inline int factoredHasRow(uint32_t ref, int row4) {
   uint32_t lowPair = ref & 0xFFFF, highPair = ref >> 16;
   return (lowHalfMasks[(lowPair << lowResultBits) | (uint32_t)(row4 & ((1 << lowResultBits) - 1))]
         & highHalfMasks[(highPair << highResultBits) | (uint32_t)(row4 >> lowResultBits)]) != 0;
}
static inline const uint16_t *factoredRows(uint32_t ref, int row4, int *n, uint16_t *buf) {
   uint32_t lowPair = ref & 0xFFFF, highPair = ref >> 16;
   const int row4Low = row4 & ((1 << lowResultBits) - 1), row4High = row4 >> lowResultBits;
   int overlapMask = lowHalfMasks[(lowPair << lowResultBits) | row4Low] & highHalfMasks[(highPair << highResultBits) | row4High];
   int count = 0;
   if (overlapMask) {
      const uint16_t *lowOffsets = lowHalfOffsets + (size_t)lowPair * LOW_OFFSETS_PER_PAIR + (row4Low << 2);
      const uint16_t *highOffsets = highHalfOffsets + (size_t)highPair * HIGH_OFFSETS_PER_PAIR + (row4High << 2);
      const uint8_t *lowRows = lowHalfRows + ((size_t)lowPair << lowHalfBits);
      const uint8_t *highRows = highHalfRows + ((size_t)highPair << highHalfBits);
      const int shift = highHalfShift;
      for (int overlap = 0; overlap < 4; overlap++) {
         if (!((overlapMask >> overlap) & 1)) continue;
         for (int i = lowOffsets[overlap]; i < lowOffsets[overlap + 1]; i++) {
            const int lowRow = lowRows[i];
            for (int j = highOffsets[overlap]; j < highOffsets[overlap + 1]; j++)
               buf[count++] = (uint16_t)(lowRow | (highRows[j] << shift));
         }
      }
#ifdef FACTORED_SORTED
      for (int i = 1; i < count; i++) {      /* insertion sort into per-pair-table order */
         uint16_t x = buf[i]; uint32_t rankX = pairTableRank[x]; int j = i - 1;
         while (j >= 0 && pairTableRank[buf[j]] > rankX) { buf[j + 1] = buf[j]; j--; }
         buf[j + 1] = x;
      }
#endif
   }
   *n = count;
   return buf;
}

/* build one half-table: for every (row1 half, row2 half) a block of row3 halves */
static void buildHalfTable(int isHigh, int halfBits, int resultBits, uint8_t **rowsOut, uint16_t **offsetsOut, uint8_t **masksOut) {
   const long long numPairs = 1LL << (2 * halfBits);
   const int offsetsPerPair = (4 << resultBits) + 1;
   *rowsOut = (uint8_t*)malloc((size_t)(numPairs << halfBits));
   *offsetsOut = (uint16_t*)malloc((size_t)numPairs * offsetsPerPair * sizeof(uint16_t));
   *masksOut = (uint8_t*)calloc((size_t)(numPairs << resultBits), 1);
   if (!*rowsOut || !*offsetsOut || !*masksOut) { printf("Unable to allocate factored tables\n"); exit(1); }
   memusage += (numPairs << halfBits) + numPairs * offsetsPerPair * (long long)sizeof(uint16_t) + (numPairs << resultBits);
   int *result = (int*)malloc(sizeof(int) << halfBits);
   int *groupCount = (int*)malloc(sizeof(int) * (4 << resultBits));
   const int halfMask = (1 << halfBits) - 1;
   for (long long pair = 0; pair < numPairs; pair++) {
      const int r1 = (int)(pair >> halfBits), r2 = (int)(pair & halfMask);
      memset(groupCount, 0, sizeof(int) * (4 << resultBits));
      for (int r3 = 0; r3 <= halfMask; r3++) {
         int r4;
         if (!isHigh) r4 = evolveRowLow(r1, r2, r3, lowResultBits);
         else {
            r4 = evolveRowHigh(r1 << highHalfShift, r2 << highHalfShift, r3 << highHalfShift, highResultBits);
            if (r4 >= 0) r4 >>= lowResultBits;
         }
         result[r3] = r4;
         if (r4 >= 0) groupCount[(r4 << 2) | (isHigh ? (r3 & 3) : (r3 >> highHalfShift))]++;
      }
      uint16_t *pairOffsets = *offsetsOut + (size_t)pair * offsetsPerPair;
      int pos = 0;
      for (int k = 0; k < (4 << resultBits); k++) { pairOffsets[k] = (uint16_t)pos; pos += groupCount[k]; groupCount[k] = pairOffsets[k]; }
      pairOffsets[4 << resultBits] = (uint16_t)pos;
      uint8_t *pairRows = *rowsOut + ((size_t)pair << halfBits);
      uint8_t *pairMasks = *masksOut + ((size_t)pair << resultBits);
      for (int r3 = 0; r3 <= halfMask; r3++) {
         int r4 = result[r3];
         if (r4 < 0) continue;
         int overlap = isHigh ? (r3 & 3) : (r3 >> highHalfShift);
         pairRows[groupCount[(r4 << 2) | overlap]++] = (uint8_t)r3;
         pairMasks[r4] |= (uint8_t)(1 << overlap);
      }
   }
   free(result); free(groupCount);
}

static void buildFactoredTables(void) {
   lowHalfBits = (width >> 1) + 1; lowResultBits = lowHalfBits - 1; highHalfShift = lowHalfBits - 2; highHalfBits = width - highHalfShift; highResultBits = width - lowResultBits;
   buildHalfTable(0, lowHalfBits, lowResultBits, &lowHalfRows, &lowHalfOffsets, &lowHalfMasks);
   buildHalfTable(1, highHalfBits, highResultBits, &highHalfRows, &highHalfOffsets, &highHalfMasks);
   pairTableRank = (uint32_t*)malloc(sizeof(uint32_t) << width);
   for (int i = 0; i < 1 << width; i++) pairTableRank[rowOrder[i]] = (uint32_t)i;
}

/* decide the representation for this search (called once width is known) */
#ifndef SPLIT_AUTO_WIDTH
#define SPLIT_AUTO_WIDTH 12   /* auto: split lookahead from this width (measured crossover) */
#endif
void chooseTableMode(void) {
   useFactored = (tableMode == 2) || (tableMode == 0 && width >= 10);
   if (useFactored && width < 4) {
      if (tableMode == 2) fprintf(stderr, "Warning: factored tables need width >= 4; using per-pair tables.\n");
      useFactored = 0;
   }
   /* lookahead method: auto = split from SPLIT_AUTO_WIDTH, cone below; split and hybrid
   ** need the half-tables (width >= 4), which are tiny and built whenever needed */
   if (lookaheadMethod == 1) useSplitLookahead = 0;
   else if (lookaheadMethod == 2) useSplitLookahead = 1;
   else if (lookaheadMethod == 3) useSplitLookahead = 2;
   else useSplitLookahead = (width >= SPLIT_AUTO_WIDTH) ? 1 : 0;
   if (useSplitLookahead && width < 4) useSplitLookahead = 0;
}

#define PAIR_TABLE_REF(r1, r2)            ((TableRef)tableFor((r1), (r2)))
#define PAIR_HAS_ROW(t, r)              BITSET((const uint16_t *)(t), (r))
#define PAIR_ROWS(t, r, n, b)        bucketEntries((const uint16_t *)(t), (r), (n), (b))
#define FACTORED_TABLE_REF(r1, r2)            ((TableRef)factoredTableRef((r1), (r2)))
#define FACTORED_HAS_ROW(t, r)              factoredHasRow((uint32_t)(t), (r))
#define FACTORED_ROWS(t, r, n, b)        factoredRows((uint32_t)(t), (r), (n), (b))
#define TABLE_REF(fac, r1, r2)       ((fac) ? FACTORED_TABLE_REF(r1, r2) : PAIR_TABLE_REF(r1, r2))
#define TABLE_HAS_ROW(fac, t, r)         ((fac) ? FACTORED_HAS_ROW(t, r) : PAIR_HAS_ROW(t, r))
#define TABLE_ROWS(fac, t, r, n, b)   ((fac) ? FACTORED_ROWS(t, r, n, b) : PAIR_ROWS(t, r, n, b))

/* Given rows row1, row2, and row3, getSuccessors() gives the list (p) *\
** of all rows XXXX such that                                          **
**                                                                     **
**    row1   evolves into                                              **
**    row2 ----------------> row3                                      **
**    XXXX                                                             **
**                                                                     **
\* and its length (n).  buf must hold DECODE_CAP entries.              */
void getSuccessors(int row1, int row2, int row3, const uint16_t** p, int *n, uint16_t *buf) {
   if (useFactored) *p = FACTORED_ROWS(FACTORED_TABLE_REF(row1, row2), row3, n, buf);
   else             *p = PAIR_ROWS(PAIR_TABLE_REF(row1, row2), row3, n, buf);
}

/* Same, and also returns a 64-bit identity of the list (for the lookahead cache). */
void getSuccessorsWithId(int row1, int row2, int row3, const uint16_t** p, int *n, uint16_t *buf, uint64_t *id) {
   if (useFactored) {
      uint32_t t = factoredTableRef(row1, row2);
      *id = (((uint64_t)t << 16) | (uint64_t)row3) + 1;   /* non-zero, unique per list */
      *p = factoredRows(t, row3, n, buf);
      return;
   }
   const uint16_t *t = tableFor(row1, row2);
   int start;
   *n = bucketFind(t, row3, &start);
   *id = ((uint64_t)(uintptr_t)t << 16) | (uint64_t)start;
   *p = *n ? bucketRead(t, start, *n, buf) : buf;
}

uint8_t *causesBirth;
row *flipRow;
int *buildScratch;       /* buildScratch to be parceled out between threads */
/* tableDedupe: open-addressing set (linear probing) of the distinct lookup tables  */
/* built so far; entry = row12 of a representative table plus one, 0 = empty.    */
/* It is sized by the number of tables actually built, not by 2^(2*width).       */
int *tableDedupe;
size_t tableDedupeCap = 0, tableDedupeCount = 0;
uint16_t *tableScratch;    /* per-thread scratch for building one table before it is stored */
uint16_t *rowOrder;
void genStatCounts(void);

void makeFlip(void) {
	row theRow;
	int i;
	for (theRow = 0; theRow < (1<<width); theRow++) {
		row flippedRow = 0;
		for (i = 0; i < width; i++)
			if (theRow & (1<<i))
				flippedRow |= 1 << (width - i - 1);
		flipRow[theRow] = flippedRow;
	}
}

void makeTables(void) {
   flipRow = (row*)malloc(sizeof(*flipRow)<<width);
   makeFlip();
   causesBirth = (uint8_t*)calloc(1LL<<width, sizeof(*causesBirth));
   if (!useFactored) {
      pairTable = (_Atomic(uint16_t*) *)calloc(1LL<<(width*2), sizeof(*pairTable));
      /* pairTable comes from calloc(), so every entry is already NULL.  It is deliberately */
      /* not initialized by a loop: the untouched pages are never committed to RAM.   */
      tableDedupeCap = 1024;
      tableDedupe = (int *)calloc(tableDedupeCap, sizeof(*tableDedupe));
      tableDedupeCount = 0;
      memusage += (sizeof(*pairTable) << (width*2)) + tableDedupeCap * sizeof(*tableDedupe);
   }
   rowOrderWeight = (uint32_t *)calloc(1LL << width, sizeof(*rowOrderWeight));
   
   for (uint32_t j = 0; j < 1LLU << width; ++j)
      causesBirth[j] = (evolveRow(j,0,0) ? 1 : 0);
   for (uint32_t j = 0; j < 1LLU << width; ++j)
      rowOrderWeight[j] = 0;
   if (!useFactored) {
      buildScratch = (int *)calloc((4LL*params[P_NUMTHREADS])<<width, sizeof(int));
      tableScratch = (uint16_t *)calloc((long long)params[P_NUMTHREADS] * SCRATCH_WORDS, sizeof(uint16_t));
   }
   if (params[P_REORDER] == 1)
      genStatCounts();
   if (params[P_REORDER] == 2)   /* this option currently cannot be set at runtime */
      for (int i=1; i<1<<width; i++)
         rowOrderWeight[i] = 1 + rowOrderWeight[i & (i - 1)];
   rowOrderWeight[0] = 0xffffffff;    /* Maximum value so empty row is chosen first */
   rowOrder = (uint16_t *)calloc(1LL << width, sizeof(uint16_t));
   for (int i=0; i<1<<width; i++)
      rowOrder[i] = (uint16_t) ((1<<width)-1-i);
   if (params[P_REORDER] != 0)
      sortRows(rowOrder, 1<<width);
   if (useFactored || useSplitLookahead)
      buildFactoredTables();
   if (!useFactored)
      for (int row2=0; row2<1<<width; row2++)
         buildTable(0, row2);
}

uint16_t *tableChunk;
long long int tableChunkLeft = 0;

/* reduce fragmentation by allocating chunks larger */
/* than needed and parceling out the small pieces.  */
/* Chunks are large (at least 8 MB) and, on Linux, 2 MB aligned with transparent  */
/* huge pages requested: table construction otherwise spends a large share of its */
/* time in page faults (one per 4 KB of table).                                   */
#define TABLE_CHUNK_MIN_WORDS ((long long)((8LL<<20) / sizeof(uint16_t)))
uint16_t *tableAlloc(int siz) {
   if (siz + (1<<width) > tableChunkLeft) {
      long long need = (1LL << (2 * width)) + (1LL<<width) + MASKWORDS + PREWORDS + 8;
      tableChunkLeft = need < TABLE_CHUNK_MIN_WORDS ? TABLE_CHUNK_MIN_WORDS : need;
#ifdef __linux__
      {
         const size_t HUGE = (size_t)2 << 20;
         size_t bytes = ((size_t)tableChunkLeft * sizeof(uint16_t) + HUGE - 1) & ~(HUGE - 1);
         tableChunkLeft = (long long)(bytes / sizeof(uint16_t));
         void *raw = mmap(NULL, bytes + HUGE, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
         if (raw == MAP_FAILED) {
            tableChunk = (uint16_t *)calloc((size_t) tableChunkLeft, sizeof(uint16_t));
         } else {
            uintptr_t al = ((uintptr_t)raw + HUGE - 1) & ~(uintptr_t)(HUGE - 1);
            tableChunk = (uint16_t *)al;
            #ifdef MADV_HUGEPAGE
            madvise(tableChunk, bytes, MADV_HUGEPAGE);
            #endif
         }
      }
#else
      tableChunk = (uint16_t *)calloc((size_t) tableChunkLeft, sizeof(uint16_t));
#endif
      if (!tableChunk) {
         printf("Unable to allocate lookup table memory\n");
         exit(1);
      }
      memusage += 2*tableChunkLeft;
      if (params[P_MEMLIMIT] >= 0 && memusage > memlimit) {
         printf("Aborting due to excessive memory usage\n");
         exit(1);
      }
   }
   uint16_t *r = tableChunk;
   tableChunk += siz;
   tableChunkLeft -= siz;
   return r;
}


/* hashTableWords() is used to identify if we've already built an identical part of */
/* the lookup table.  Four independent 64-bit lanes keep the multiplier busy   */
/* (the previous byte-serial h = h*3 + x was a latency-bound chain).  Only the */
/* distribution matters here, never the value, so the table contents and the   */
/* search are unaffected by the choice of hash function.                       */
static inline uint32_t hashTableWords(const uint16_t *theRow, int siz) {
   const unsigned char *p = (const unsigned char *)theRow;
   size_t n = (size_t)siz * sizeof(uint16_t);
   const uint64_t K = 0x9E3779B97F4A7C15ULL;
   uint64_t h0 = 0x243F6A8885A308D3ULL ^ n, h1 = 0x13198A2E03707344ULL, h2 = 0xA4093822299F31D0ULL, h3 = 0x082EFA98EC4E6C89ULL;
   size_t i = 0;
   for (; i + 32 <= n; i += 32) {
      uint64_t v0, v1, v2, v3;
      memcpy(&v0, p + i, 8); memcpy(&v1, p + i + 8, 8);
      memcpy(&v2, p + i + 16, 8); memcpy(&v3, p + i + 24, 8);
      h0 = (h0 ^ v0) * K; h0 ^= h0 >> 29;
      h1 = (h1 ^ v1) * K; h1 ^= h1 >> 29;
      h2 = (h2 ^ v2) * K; h2 ^= h2 >> 29;
      h3 = (h3 ^ v3) * K; h3 ^= h3 >> 29;
   }
   for (; i < n; i += 8) {
      uint64_t v = 0;
      memcpy(&v, p + i, (n - i) < 8 ? (n - i) : 8);
      h0 = (h0 ^ v) * K; h0 ^= h0 >> 29;
   }
   uint64_t h = (h0 ^ (h1 * 3) ^ (h2 * 5) ^ (h3 * 7)) * K;
   h ^= h >> 32;
   return (uint32_t)h;
}

/* Called with the updateTable lock held: double tableDedupe and reinsert every */
/* representative table.                                                    */
static void growTableDedupe(void) {
   size_t ncap = tableDedupeCap * 2;
   int *nh = (int *)calloc(ncap, sizeof(*nh));
   for (size_t i = 0; i < tableDedupeCap; i++) {
      int v = tableDedupe[i];
      if (!v) continue;
      const uint16_t *t = pairTable[v - 1];
      size_t h = hashTableWords(t, tableWords(t)) & (ncap - 1);
      while (nh[h]) h = (h + 1) & (ncap - 1);
      nh[h] = v;
   }
   free(tableDedupe);
   memusage += (long long)(ncap - tableDedupeCap) * sizeof(*nh);
   tableDedupe = nh;
   tableDedupeCap = ncap;
}

uint16_t *buildTable(int row1, int row2) {
   const int N = 1 << width;
   const int maskWords = MASKWORDS;
   const int tid = omp_get_thread_num();
   /* Per-thread scratch: row4List = row4 list, row3List = matching row3 list,       */
   /* bucketCount = per-row4 counts (later running write positions), partialRow4 = partial   */
   /* evolution tables.                                                         */
   int *row4List  = buildScratch + ((4LL * tid) << width);
   int *row3List = row4List + N;
   int *bucketCount    = row3List + N;
   int *partialRow4 = bucketCount + N;
   uint16_t *tbl = tableScratch + (long long)tid * SCRATCH_WORDS;
   int numEntries = 0;
   
   /* For each row3, find all row4 such that   *\
   **                                          **
   **      row1   evolves into                 **
   **      row2 ----------------> row4         **
   **      row3                                **
   **                                          **
   ** list of row4s is stored in row4List and     **
   \* corresponding row3s are stored in row3List */
   memset(bucketCount, 0, (size_t)N * sizeof(int));
   if (width < 4) {
      for (int row3=0; row3<N; row3++)
         partialRow4[row3] = evolveRow(row1, row2, row3);
      for (int row3i = 0; row3i < N; row3i++) {
         int row3 = rowOrder[row3i];
         int row4 = partialRow4[row3];
         if (row4 < 0)
            continue;
         row3List[numEntries] = row3;
         row4List[numEntries++] = row4;
         bucketCount[row4]++;
      }
   } else {
      int lowbitcount = (width >> 1) + 1;
      int hibitcount = ((width + 1) >> 1) + 1;
      int hishift = lowbitcount - 2;
      int lowcount = 1 << lowbitcount;
      int *lowTab = partialRow4;
      int *hiTab = partialRow4 + lowcount;
      for (int row3=0; row3<lowcount; row3++)
         lowTab[row3] = evolveRowLow(row1, row2, row3, lowbitcount-1);
      for (int row3=0; row3<N; row3 += 1<<hishift)
         hiTab[row3>>hishift] = evolveRowHigh(row1, row2, row3, hibitcount-1);
      /* row4 = low part | high part; either part is -1 if the row is impossible */
      const int lm = lowcount - 1;
      for (int row3i = 0; row3i < N; row3i++) {
         int row3 = rowOrder[row3i];
         int row4 = lowTab[row3 & lm] | hiTab[row3 >> hishift];
         if (row4 < 0)
            continue;
         row3List[numEntries] = row3;
         row4List[numEntries++] = row4;
         bucketCount[row4]++;
      }
   }
   
   /* Build the table in a private buffer (layout described at tabEntries()).    */
   const int prefixWords = PREWORDS;
   int numBuckets = 0;
   for (int row4 = 0; row4 < N; row4++) numBuckets += (bucketCount[row4] != 0);
   const int offBase = maskWords + prefixWords + 2;
   const int dataBase = offBase + numBuckets + 1;
#ifdef PACKED
   const int dataWords = ((numEntries * width + 15) >> 4) + 8;
   const int tableSize = dataBase + dataWords;
   const int firstPos = 0;
#else
   const int tableSize = dataBase + numEntries;
   const int firstPos = dataBase;
#endif
   {
      int start = firstPos, idx = 0;
      for (int j = 0; j < prefixWords; j++) {
         uint64_t m = 0;
         tbl[maskWords + j] = (uint16_t)idx;
         const int lim = (N - (j << 6) < 64) ? N - (j << 6) : 64;
         for (int bit = 0; bit < lim; bit++) {
            int c = bucketCount[(j << 6) + bit];
            if (c) {
               m |= 1ULL << bit;
               tbl[offBase + idx++] = (uint16_t)start;
               bucketCount[(j << 6) + bit] = start;        /* now the next write position */
               start += c;
            }
         }
         for (int q = 0; q < 4; q++) tbl[(j << 2) + q] = (uint16_t)(m >> (16 * q));
      }
      tbl[offBase + numBuckets] = (uint16_t)start;
      tbl[maskWords + prefixWords] = (uint16_t)tableSize;
      tbl[maskWords + prefixWords + 1] = (uint16_t)dataBase;
   }
#ifdef PACKED
   memset(tbl + dataBase, 0, (size_t)dataWords * sizeof(uint16_t));
   for (int i = 0; i < numEntries; i++) {
      unsigned bp = (unsigned)(bucketCount[row4List[i]]++) * width;
      uint32_t v = (uint32_t)row3List[i] << (bp & 15);
      tbl[dataBase + (bp >> 4)] |= (uint16_t)v;
      tbl[dataBase + (bp >> 4) + 1] |= (uint16_t)(v >> 16);
   }
#else
   for (int i = 0; i < numEntries; i++)
      tbl[bucketCount[row4List[i]]++] = (uint16_t)row3List[i];
#endif
   
   const uint32_t hash32 = hashTableWords(tbl, tableSize);
   
   /* tableAlloc, and all operations that read from or write to the stored tables, */
   /* tableDedupe, and pairTable must be included in a critical region.                 */
   uint16_t *result;
   #pragma omp critical(updateTable)
   {
      size_t h = hash32 & (tableDedupeCap - 1);
      for (;;) {
         int v = tableDedupe[h];
         if (v == 0) {
            /* a new distinct table: store it and register it */
            uint16_t *tableStart = tableAlloc(tableSize);
            memcpy(tableStart, tbl, (size_t)tableSize * sizeof(uint16_t));
            result = tableStart;
            tableDedupe[h] = (row1 << width) + row2 + 1;
            if (++tableDedupeCount * 2 > tableDedupeCap) {
               /* pairTable must already point at it before it can be rehashed */
               atomic_store_explicit(&pairTable[(row1<<width)+row2], result, memory_order_release);
               growTableDedupe();
            }
            break;
         }
         /* Maybe two different row12s result in the exact same rows for the */
         /* lookup table. This prevents two different threads from trying to */
         /* build the same part of the lookup table.                         */
         const uint16_t *other = atomic_load_explicit(&pairTable[v - 1], memory_order_relaxed);
         if (tableWords(other) == tableSize &&
             memcmp(tbl, other, (size_t)tableSize * sizeof(uint16_t)) == 0) {
            result = (uint16_t *)other;
            break;
         }
         h = (h + 1) & (tableDedupeCap - 1);
      }
      atomic_store_explicit(&pairTable[(row1<<width)+row2], result, memory_order_release);
   }
   
   return result;
}

/*   We calculate the stats using a 2 * 64 << width array.  We use a
**   leading 1 to separate them.  Index 1 aaa bb cc dd represents
**   the count for a result of aaa when the last two bits of row1, row2,
**   and row3 were bb, cc, and dd, respectively.  We have to manage
**   the edge conditions appropriately.
*/
void genStatCounts(void) {
   /* unsigned: these counts overflow 32 bits at width 14 (signed overflow was undefined */
   /* behaviour); wrapping gives the same values the original computed in practice.    */
   uint32_t *cnt = (uint32_t*)calloc(1LL << width, 128 * sizeof(uint32_t));
   for (int i=0; i<128<<width; i++)
      cnt[i] = 0;
   int s = 0;
   if (params[P_SYMMETRY] == SYM_ODD)
      s = 2;
   else if (params[P_SYMMETRY] == SYM_EVEN)
      s = 1;
   else
      s = width + 2;
   /* left side: never permit generation left of row4 */
   for (int row1=0; row1<2; row1++)
      for (int row2=0; row2<2; row2++)
         for (int row3=0; row3<2; row3++)
            if (evolveBit(row1, row2, row3) == 0)
               cnt[(1<<6) + (row1 << 4) + (row2 << 2) + row3]++;
   for (int nb=0; nb<width; nb++) {
      for (int row1=0; row1<8; row1++)
         for (int row2=0; row2<8; row2++)
            for (int row3=0; row3<8; row3++) {
               if (nb == width-1)
                  if ((((row1 >> s) ^ row1) & 1) ||
                      (((row2 >> s) ^ row2) & 1) ||
                      (((row3 >> s) ^ row3) & 1))
                     continue;
               int row4b = evolveBit(row1, row2, row3);
               for (int row4=0; row4<1<<nb; row4++)
                  cnt[(((((1<<nb) + row4) << 1) + row4b) << 6) +
                    ((row1 & 3) << 4) + ((row2 & 3) << 2) + (row3 & 3)] +=
                     cnt[(((1<<nb) + row4) << 6) +
                       ((row1 >> 1) << 4) + ((row2 >> 1) << 2) + (row3 >> 1)];
            }
   }
   /* right side; check left, and accumulate into rowOrderWeight */
   for (int row1=0; row1<4; row1++)
      for (int row2=0; row2<4; row2++)
         for (int row3=0; row3<4; row3++)
            if (params[P_SYMMETRY] != SYM_ASYM ||
                evolveBit(row1<<1, row2<<1, row3<<1) == 0)
               for (int row4=0; row4<1<<width; row4++)
                  rowOrderWeight[row4] +=
                     cnt[(((1<<width) + row4) << 6) +
                       (row1 << 4) + (row2 << 2) + row3];
   free(cnt);
}

/* ====================================================== */
/*  Hash table for detecting equivalent partial patterns  */
/* ====================================================== */

void resetHash(void) {
   if (seenHash != 0) memset(seenHash,0,4*HASHSIZE);
}


static inline long seenHashIndex(node b, row r) {
   /* Unsigned arithmetic (signed overflow was undefined behaviour).  The right    */
   /* shifts stay arithmetic so the hash values are exactly the ones the original */
   /* produced on 64-bit platforms.                                               */
   uint64_t h = r;
   if (params[P_SYMMETRY] == SYM_ASYM) h += flipRow[r];
   int i;
   for (i = 0; i < nRowsInState; i++) {
      h = (h * 269) + ROW(b);
      if (params[P_SYMMETRY] == SYM_ASYM) h += flipRow[ROW(b)];
      b = PARENT(b);
   }
   h += (uint64_t)((int64_t)h >> 16) * 269;
   h += (uint64_t)((int64_t)h >> 8) * 269;
   return (long)(h & HASHMASK);
}

/* test if q+r is same as p */
static inline int sameState(node p, node q, row r) {
   int i;
   for (i = 0; i < nRowsInState; i++) {
      if (p >= QSIZE || q >= QSIZE || EMPTY(p) || EMPTY(q)) return 0;   /* sanity check */
      if (ROW(p) != r) return 0;
      p = PARENT(p);
      r = ROW(q);
      q = PARENT(q);
   }
   return 1;
}

static inline int sameStateFlipped(node p, node q, row r) {
   int i;
   for (i = 0; i < nRowsInState; i++) {
      if (p >= QSIZE || q >= QSIZE || EMPTY(p) || EMPTY(q)) return 0;   /* sanity check */
      if (flipRow[ROW(p)] != r) return 0;
      p = PARENT(p);
      r = ROW(q);
      q = PARENT(q);
   }
   return 1;
}

/* test if we've seen this child before */
static inline int isVisited(node b, row r) {
   if (sameState(0,b,r)) return 1;
   if (seenHash != 0) {
      int hashVal = seenHashIndex(b,r);
      node hashNode = seenHash[hashVal];
      if (hashNode == 0)
         return 0;
      else if (sameState(hashNode,b,r))
         return 1;
      else if (params[P_SYMMETRY] == SYM_ASYM && sameStateFlipped(hashNode,b,r))
         return 1;
   }
   return 0;
}

/* Add node (NOT child) to hash table */
static inline void setVisited(node b) {
   if (seenHash != 0) seenHash[ seenHashIndex(PARENT(b),ROW(b)) ] = b;
}

/* ============================================== */
/*  Output patterns found by successful searches  */
/* ============================================== */

#define MAXRLELINEWIDTH 63
int RLEcount = 0;
int RLElineWidth = 0;
char RLEchar;
char * patternBuf;

void bufRLE(char c) {
  if (RLEcount > 0 && c != RLEchar) {
      if (RLElineWidth++ >= MAXRLELINEWIDTH) {
         if (RLEchar != '\n') sprintf(patternBuf+strlen(patternBuf),"\n");
         RLElineWidth = 0;
      }
    if (RLEcount == 1) strncat(patternBuf,&RLEchar,1);
    else {
       sprintf(patternBuf+strlen(patternBuf),"%d%c", RLEcount, RLEchar);
       RLElineWidth ++;
       if (RLEcount > 9) RLElineWidth++;
    }
    RLEcount = 0;
    if (RLEchar == '\n') RLElineWidth = 0;
  }
  if (c != '\0') {
    RLEcount++;
    RLEchar = c;
  } else RLElineWidth = 0;
}

void bufRow(unsigned long rr, unsigned long r, int shift) {
   while (r | rr) {
      if (shift == 0)
         bufRLE(r & 1 ? 'o' : 'b');
      else shift--;
      r >>= 1;
      if (rr & 1) r |= (1UL<<31);
      rr >>= 1;
   }
   bufRLE('$');
}

/* Avoid Intel shift bug */
static inline unsigned long
safeShift(unsigned long r, int i) {
    unsigned long rr = r;
    while (i>16) { rr >>= 16; i-=16;}
    return (rr>>i);
}

int outRowsAlloc =  0;
unsigned long * outRowsLowBuf;
unsigned long * outRowsHighBuf;
int prevOutRowCount = 0;
unsigned long * prevOutRowsLow;
unsigned long * prevOutRowsHigh;

/* Buffers RLE into patternBuf; returns 1 if pattern successfully buffered */ 
int bufferPattern(node b, row *pathRows, int nodeRow, uint32_t lastRow, int printExpected) {
   node c;
   int numRows = 0;
   int patternWidth;
   int outRowsNeeded;
   int p, i, j, margin;
   unsigned long *rowsLow, *rowsHigh;

   uint32_t currRow = lastRow;
   int nDeepRows = 0;
   int nodeDiff;

   if (pathRows != NULL){
      while (pathRows[currRow] == 0){
         if (currRow == 0){
            if (!printExpected) return 0;
            printf("Success called on search root!\n");
            aborting = 1;
            return 0;
         }
         currRow--;
      }
      nDeepRows = (currRow / period) - 1;
      nodeDiff = nodeRow - period - (currRow%period);
      nodeRow -= nodeDiff;

      for (j = 0; j < nodeDiff; j++){
         b = PARENT(b);
      }
      currRow = currRow - period + 1;
      numRows = nDeepRows;
      
   }
   else{
      /* shift until we find nonzero row.
         then shift PERIOD-1 more times to get leading edge of glider */
      while (ROW(b) == 0) {
         b = PARENT(b);
         if (b == 0) {
            if (!printExpected) return 0;
            printf("Success called on search root!\n");
            aborting = 1;
            return 0;
         }
      }
   }
   if (numRows < 0) numRows = 0;
   
   for (p = 0; p < period-1; p++) b = PARENT(b);
   if (b == 0) {
      if (!printExpected) return 0;
      printf("Success called on search root!\n");
      aborting = 1;
      return 0;
   }
   
   /* count rows */
   c = b;
   while (c != 0) {
      for (p = 0; p < period; p++)
         c = PARENT(c);
      numRows++;
   }
   
   /* build data structure of rows so we can reduce width etc */
   outRowsNeeded = numRows+MAXWIDTH+1;
   if (!outRowsAlloc)
   {
      outRowsAlloc = outRowsNeeded;
      outRowsLowBuf = (unsigned long*)malloc(outRowsAlloc * sizeof(unsigned long));
      outRowsHighBuf = (unsigned long*)malloc(outRowsAlloc * sizeof(unsigned long));
      prevOutRowsLow = (unsigned long*)calloc(outRowsAlloc, sizeof(unsigned long));
      prevOutRowsHigh = (unsigned long*)calloc(outRowsAlloc, sizeof(unsigned long));
      patternBuf = (char*)malloc(((2 * MAXWIDTH + 4) * outRowsAlloc + 300) * sizeof(char));
   }
   else if (outRowsAlloc < outRowsNeeded)
   {
      outRowsAlloc = outRowsNeeded;
      outRowsLowBuf = (unsigned long*)realloc(outRowsLowBuf, outRowsAlloc * sizeof(unsigned long));
      outRowsHighBuf = (unsigned long*)realloc(outRowsHighBuf, outRowsAlloc * sizeof(unsigned long));
   }
   rowsLow  = outRowsLowBuf;
   rowsHigh = outRowsHighBuf;
   
   for (i = 0; i <= numRows+MAXWIDTH; i++) rowsLow[i]=rowsHigh[i]=0;
   for (i = numRows - 1; i >= 0; i--) {
      uint32_t r;
      if (nDeepRows > 0){
         r = pathRows[currRow];
         currRow -= period;
         nDeepRows--;
      }
      else{
         r = ROW(b);
         for (p = 0; p < period; p++) {
            b = PARENT(b);
         }
      }
      switch(params[P_SYMMETRY]) {
         case SYM_ASYM:
            rowsLow[i] = r;
            break;

         case SYM_ODD:
            rowsLow[i] = r << (MAXWIDTH - 1);
            rowsHigh[i] = r >> (32 - (MAXWIDTH - 1));
            for (j = 1; j < MAXWIDTH; j++)
               if (r & (1<<j))
                  rowsLow[i] |= 1 << (MAXWIDTH - 1 - j);
            break;

         case SYM_EVEN:
            rowsLow[i] = r << MAXWIDTH;
            rowsHigh[i] = r >> (32 - MAXWIDTH);
            for (j = 0; j < MAXWIDTH; j++)
               if (r & (1<<j))
                  rowsLow[i] |= 1 << (MAXWIDTH - 1 - j);
            break;

         case SYM_GUTTER:
            rowsLow[i] = r << (MAXWIDTH + 1);
            rowsHigh[i] = r >> (32 - (MAXWIDTH + 1));
            for (j = 0; j < MAXWIDTH; j++)
               if (r & (1<<j))
                  rowsLow[i+gutterSkew] |= 1 << (MAXWIDTH - 1 - j);
            break;

         default:
            printError("unexpected symmetry type in success()");
            return 0;
      }
   }
   
   /* normalize nrows to only include blank rows */
   numRows += MAXWIDTH;
   while (numRows>0 && rowsLow[numRows-1] == 0 && rowsHigh[numRows-1] == 0) numRows--;
   while (rowsLow[0] == 0 && rowsHigh[0] == 0 && numRows>0) {
      rowsLow++;
      rowsHigh++;
      numRows--;
   }
   
   /* sanity check: is the pattern nonempty? */
   int allEmpty = 1;
   for (i = 0; i < numRows; i++){
      if (rowsLow[i]){
         allEmpty = 0;
         break;
      }
   }
   if (allEmpty){
      free(patternBuf);
      patternBuf = NULL;
      return 0;
   }
   
   /* make at least one row have nonzero first bit */
   i = 0;
   while ((rowsLow[i] & 1) == 0) {
      for (i = 0; i < numRows && (rowsLow[i] & 1) == 0; i++) { }
      if (i == numRows) {
         for (i = 0; i < numRows; i++) {
            rowsLow[i] >>= 1;
            if (rowsHigh[i] & 1) rowsLow[i] |= (1UL<<31);
            rowsHigh[i] >>= 1;
         }
         i = 0;
      }
   }
   
   patternWidth = 0;
   for (i = 0; i < numRows; i++)
      while (safeShift(rowsHigh[i],patternWidth))
         patternWidth++;
   if (patternWidth) patternWidth += 32;
   for (i = 0; i < numRows; i++)
      while (safeShift(rowsLow[i],patternWidth))
         patternWidth++;
   
      
   /* compute margin on other side of width */
   margin = 0;

   /* make sure we didn't just output the exact same pattern */
   if (printExpected){
      if (numRows == prevOutRowCount) {
         int different = 0;
         for (i = 0; i < numRows && !different; i++)
            different = (rowsLow[i] != prevOutRowsLow[i] || rowsHigh[i] != prevOutRowsHigh[i]);
         if (!different) return 0;
      }
      
      /* replace previous saved rows with new rows */
      prevOutRowCount = numRows;
      prevOutRowsLow = (unsigned long*)realloc(prevOutRowsLow, outRowsAlloc * sizeof(unsigned long));
      prevOutRowsHigh = (unsigned long*)realloc(prevOutRowsHigh, outRowsAlloc * sizeof(unsigned long));
      memcpy(prevOutRowsLow, rowsLow, numRows * sizeof(unsigned long));
      memcpy(prevOutRowsHigh, rowsHigh, numRows * sizeof(unsigned long));
   }
   
   /* Buffer output */
   patternBuf = (char*)realloc(patternBuf, ((2 * MAXWIDTH + 4) * outRowsAlloc + 300) * sizeof(char));
   
   sprintf(patternBuf,"x = %d, y = %d, rule = %s\n", patternWidth - margin, numRows, baseRule);
   
   int theBufRow = -1;
   while (theBufRow++ < numRows){
      bufRow(rowsHigh[theBufRow], rowsLow[theBufRow], 0);
   }
   RLEcount = 1;     /* prevents erroneous printing of '2' at end of RLE */
   RLEchar = '!';
   bufRLE('\0');
   sprintf(patternBuf+strlen(patternBuf),"\n");
   
   if (printExpected){
      numFound++;
      if (params[P_NUMSHIPS] > 0){
         if (--params[P_NUMSHIPS] == 0) aborting = 3;  /* use 3 to flag that we reached ship limit */
      }
   }
   
   return 1;
}

#ifndef QSIMPLE
void makeSubperiodTables(void);
int subperiodic(node x, row *pathRows, int nodeRow, uint32_t lastRow);
#endif

void success(node b, row *pathRows, int nodeRow, uint32_t lastRow) {
#ifndef QSIMPLE
   if (subperiodic(b, pathRows, nodeRow, lastRow)) return;
#endif
   if (bufferPattern(b, pathRows, nodeRow, lastRow, 1))
      printf("\n%s\n",patternBuf);
   fflush(stdout);
}

/* Test if this is a node at which we can stop */
int terminal(node n) {
   int p;

   for (p = 0; p < period; p++) {   /* last row in each phase must be zero */
      if (ROW(n) != 0) return 0;
      n = PARENT(n);
   }

   for (p = 0; p < period; p++) {
      if (causesBirth[ROW(n)]) return 0;
      n = PARENT(n);
   }
   return 1;
}

/* ================================================ */
/*  Queue of partial patterns still to be examined  */
/* ================================================ */

node qHead,qTail;

/* Queue dimensions required during save/restore */
node qStart; /* index of first node in queue */
node qEnd;   /* index of first unused node after end of queue */

/* Maintain the phase of queue nodes (which of the `period` generations a node's last
** row belongs to).  All nodes of one BFS generation share a phase; queuePhase is the
** phase of the generation being dequeued and nextRephase is the first node of the
** next generation, at which dequeue() advances queuePhase.  After the queue is
** compacted this must be reinitialized by rephase().  peekPhase() gives the phase of
** any node still in the queue.  (A global copy for the last dequeued node, `phase`,
** was never read and has been removed.)
*/
int queuePhase = 0;
node nextRephase = 0;

void rephase(void) {
   node x, y;
   while (qHead < qTail && EMPTY(qHead)) qHead++;   /* skip empty queue cells */
   x = qHead;   /* find next item in queue */
   queuePhase = period - 1;
   while (x != 0) {
      x = PARENT(x);
      queuePhase++;
   }
   queuePhase %= period;

   /* Now walk forward through queue finding breakpoints between each generation. */
   /* invariants: y is always the first in its generation                         */
   x = 0; y = 0;
   while (y <= qHead) {
      ++x;
      if (x >= qTail || (!EMPTY(x) && PARENT(x) >= y)) y = x;
   }
   nextRephase = y;
}

/* peekPhase() returns the phase of an element in the queue.  This only */
/* works for queue elements and should NOT be used on other nodes.      */
int peekPhase(node i) {
   return (i < nextRephase? queuePhase : (queuePhase+1)%period);
}

/* Test queue status */
static inline int qIsEmpty(void) {
   while (qHead < qTail && EMPTY(qHead)){
      ++qHead;
      ++extQHead;
   }
   return (qTail == qHead);
}

void qFull(void) {
    if (aborting != 2) {
      printf("Exceeded %llu node limit, search aborted\n", QSIZE);
      fflush(stdout);
      aborting = 2;
   }
}

static inline void enqueue(node b, row r) {
   /* If the queue is full the search is aborted (qFull); never write past the end */
   /* of rows[] or base[].                                                          */
   if (qTail >= QSIZE) { qFull(); return; }
   node tempQTail = qTail;
   node i = qTail++;
   if (FIRSTBASE(i)) {
      queueBase[i>>BASEBITS] = b;
      queueRows[i] = r;
   } else {
      long o = b - queueBase[i>>BASEBITS];
      if (o < 0 || o >(long) MAX_PARENT_OFFSET) {   /* offset out of range */
         while (!FIRSTBASE(i)) {
            queueRows[i] = (row) -1;
            if (qTail >= QSIZE) break;
            i = qTail++;
         }
         if (FIRSTBASE(i) && i < QSIZE) {
            queueBase[i>>BASEBITS] = b;
            queueRows[i] = r;
         }
         else qFull();
      } else queueRows[i] = (row)(o << width) + r;
   }
   
   /* update tail of parallel queue, but don't set value */
   extQTail += qTail - tempQTail;
   extensionIndex[extQTail] = 0;
}

static inline node dequeue(void) {
   prevExtQHead = extQHead;  /* Save old parallel queue head for use in process() */
   while (qHead < qTail && EMPTY(qHead)){
      ++qHead;
      ++extQHead;
   }
   if (qHead >= nextRephase) {
      queuePhase = (queuePhase+1)%period;
      nextRephase = qTail;
   }
   ++extQHead;
   return qHead++;
}

/* Not used, but could be useful in debugging */
static inline void pop(void) {
   qTail--;
   while (qTail > qHead && EMPTY(qTail-1)) qTail--;
}

void resetQ(void) {
   qHead = qTail = 0; extQHead = extQTail = 0;
}

//static inline int qTop() { return qTail - 1; }

/* =================== */
/*  Dump search state  */
/* =================== */

int splitNum = 0;
int dumpNum = 1;
char dumpFile[256];
const char *dumpRoot = "dump-@time-";
char loadDumpRoot[251];    /* used for loading dump root from file */
char trueDumpRoot[251];
time_t lastDumpTime;

int dumpFlag = 0;    /* Dump status flags, possible values follow */
#define DUMPRESET   (0)
#define DUMPPENDING (1)
#define DUMPFAILURE (2)
#define DUMPSUCCESS (3)

int dumpMode;  /* separate from params[P_DUMPMODE] for splitting */
#define D_DISABLED   (0)
#define D_OVERWRITE  (1)
#define D_SEQUENTIAL (2)

void parseDumpRoot(void) {
   /* All copies are bounded: dumpRoot may come from a loaded file (up to 250 chars). */
   char tempStr[256] = {'\0'};
   char tempRule[151] = {'\0'};
   char *r, *t;
   const char *s;
   memset(trueDumpRoot, '\0', sizeof trueDumpRoot);
   
   /* replace first "@time" with hex timestamp */
   if ((s = strstr(dumpRoot, "@time"))){
      snprintf(tempStr, sizeof tempStr, "%.*s", (int)(s - dumpRoot), dumpRoot);
      if (snprintf(trueDumpRoot, sizeof trueDumpRoot, "%s%06llx%s", tempStr, (long long) time(NULL) & 0xffffff, s+5) >= (int)sizeof trueDumpRoot)
         fprintf(stderr, "Warning: dump root truncated\n");
   }
   else
      snprintf(trueDumpRoot, sizeof trueDumpRoot, "%s", dumpRoot);
   
   /* replace first "@rule" with rule string ('/' replaced with '_') */
   snprintf(tempStr, sizeof tempStr, "%s", trueDumpRoot);
   if ((t = strstr(tempStr, "@rule"))){
      snprintf(tempRule, sizeof tempRule, "%s", rule);
      r = strchr(tempRule, '/');
      if (r) *r = '_';
      *t = '\0';
      if (snprintf(trueDumpRoot, sizeof trueDumpRoot, "%s%s%s", tempStr, tempRule, t+5) >= (int)sizeof trueDumpRoot)
         fprintf(stderr, "Warning: dump root truncated\n");
   }
   
   /* replace additional occurrences of '@' with '_' */
   for (r = trueDumpRoot; *r != '\0'; ++r)
      if (*r == '@')
         *r = '_';
   
   dumpRoot = trueDumpRoot;
}

FILE * openDumpFile(void) {
   FILE * fp;
   
   if (dumpMode == D_OVERWRITE) {
      sprintf(dumpFile, "%s%s", dumpRoot, (++dumpNum)%2 ? "gold" : "blue");
      return fopen(dumpFile, "w");
   }
   else if (dumpMode == D_SEQUENTIAL){
      while (dumpNum < DUMPLIMIT) {
         sprintf(dumpFile, "%s%05d", dumpRoot, dumpNum++);
         if ((fp = fopen(dumpFile, "r")))
            fclose(fp);
         else
            return fopen(dumpFile, "w");
      }
      if (dumpNum == DUMPLIMIT){
         dumpMode = D_OVERWRITE;
         return openDumpFile();
      }
   }
   return (FILE *) 0;
}

void dumpState(void) {
   FILE * fp;
   unsigned long long i,j;
   dumpFlag = DUMPFAILURE;
   if (!(fp = openDumpFile())) return;
   fprintf(fp,"%lu\n",FILEVERSION);
   fprintf(fp,"%s\n",rule);
   if (splitNum > 1)    /* append dumpNum to new dumpRoot when splitting search state */
      fprintf(fp,"%s%05d-\n",dumpRoot,dumpNum-1);
   else
      fprintf(fp,"%s\n",dumpRoot);
   for (j = 0; j < NUM_PARAMS; ++j)
      fprintf(fp,"%d\n",params[j]);
   fprintf(fp,"%d\n",width);
   fprintf(fp,"%d\n",period);
   fprintf(fp,"%d\n",offset);
   fprintf(fp,"%d\n",lastDeep);
   if (params[P_DUMPMODE] == D_SEQUENTIAL)
      fprintf(fp,"1\n");
   else
      fprintf(fp,"%d\n",dumpNum%2);
   fprintf(fp,"%"PRIu32"\n",qHead-qStart);
   fprintf(fp,"%"PRIu32"\n",qEnd-qStart);
   for (i = qStart; i < qEnd; ++i)
      fprintf(fp,"%"PRIu16"\n",queueRows[i]);
   for (i = 0; i < QSIZE; ++i){
      if (extensionIndex[i]){
         if (extensionIndex[i] > 1){
            for (j = 0; j < extensionRows[extensionIndex[i]][0] + 1LU + 2LU; ++j){
               fprintf(fp,"%"PRIu16"\n",extensionRows[extensionIndex[i]][j]);
            }
         }
         else {
            fprintf(fp,"0\n");
            j = 0;
            while (i < QSIZE && extensionIndex[i] <= 1){
               if (extensionIndex[i] == 1) ++j;
               ++i;
            }
            fprintf(fp,"%llu\n",j);
            if (i == QSIZE) break;
            --i;
         }
      }
   }
   fclose(fp);
   dumpFlag = DUMPSUCCESS;
}

/* ================================= */
/*  Compaction of nearly full queue  */
/* ================================= */

void printCompactNumber(long unsigned n) {
   char suffix;
   if (n >= 1000000) {
      n /= 100000;
      suffix = 'M';
   } else if (n >= 1000) {
      n /= 100;
      suffix = 'k';
   } else {
      printf("%lu", n);
      return;
   }

   if (n >= 100) printf("%lu", n/10);
   else printf("%lu.%lu", n/10, n%10);
   putchar(suffix);
}

long currentDepth(void) {
   long i;
   node x;
   x = qTail - 1;
   i = 1;
   while (x != 0) {
      x = PARENT(x);
      i++;
   }
   return i;
}

/* doCompact() has two parts.  The first part compresses the
** queue.  The second part consists of the last loop which
** converts parent bits back to parent pointers.  The search
** state may be saved in between.
*/
void doCompactPart1(void) {
   node x,y;
   qEnd = qTail;
   
   /* make a pass backwards from the end finding unused nodes at or before qHead */
   x = qTail - 1;
   y = qHead - 1;
   while (y > 0) {
      /* invariants: everything after y is still active.
                     everything after x points to something after y.
                     x is nonempty and points to y or something before y.
                     so, if x doesn't point to y, y must be unused and can be removed. */
      if (!EMPTY(y)) {
         if (y > PARENT(x)) queueRows[y] = (row) -1;
         else while (EMPTY(x) || PARENT(x) == y) x--;
      }
      y--;
   }
   
   /* make a pass forwards converting parent pointers to offset from prev parent ptr. */
   /* note that after unused nodes are eliminated, all these offsets are zero or one. */
   y = 0;
   for (x = 0; x < qTail; x++) if (!EMPTY(x)) {
      if (PARENT(x) == y) queueRows[x] = ROW(x);
      else {
         y = PARENT(x);
         queueRows[x] = (1<<width) + ROW(x);
      }
   }
   
   /* Make a pass backwards compacting gaps.
   **
   ** For most times we run this, it could be combined with the next phase, but
   ** every once in a while the process of repacking the remaining items causes them
   ** to use *MORE* space than they did before they were repacked (because of the need
   ** to leave empty space when ROFFSET gets too big) and without this phase the repacked
   ** stuff overlaps the not-yet-repacked stuff causing major badness.
   ** 
   ** For this phase, y points to the current item to be repacked, and x points
   ** to the next free place to pack an item.
   */
   x = y = qTail-1;
   for (;;) {
      if (qHead == y) qHead = x;
      if (!EMPTY(y)) {
         queueRows[x] = queueRows[y];
         x--;
      }
      if (y-- == 0) break;    /* circumlocution for while (y >= 0) because x is unsigned */
   }
   qStart = ++x;     /* mark start of queue */
}

void doCompactPart2(void) {
   node x,y;
   uint32_t i, j;
   int k;
   
   /* Make a pass forwards converting parent bits back to parent pointers.
   ** 
   ** For this phase, x points to the current item to be repacked, and y points
   ** to the parent of the previously repacked item.
   ** After the previous pass, x is initialized to first nonempty item,
   ** and all items after x are nonempty. 
   */
   qTail = 0; y = 0;
   resetHash();
   for (x = qStart; x < qEnd; x++) {
      if (PARENT_OFFSET(x)) {   /* skip forward to next parent */
         y++;
         while (EMPTY(y)) y++;
      }
      enqueue(y,ROW(x));
      //if (aborting) return;    /* why is this here? value of aborting is not changed by enqueue(). */
      if (qHead == x) qHead = qTail - 1;
      setVisited(qTail - 1);
   }
   rephase();
   
   /* Repack nonzero depth-first extension indices to the end of the queue */
   j = QSIZE - 1;
   for (i = QSIZE; i > 0; --i){
      if (extensionIndex[i - 1]){
         extensionIndex[j] = extensionIndex[i - 1];
         extensionIndex[i - 1] = 0;
         --j;
      }
   }
   
   /* Sanity check: extension queue should not take up all available space */
   if (extensionIndex[0]){
      fprintf(stderr,"Error: extension queue has too many elements.\n");
      exit(1);
   }
   
   /* Respace depth-first extension queue to match node queue */
   i = 0;
   j = 0;
   while (j < QSIZE && !extensionIndex[j]) ++j;
   for (x = qHead; x < qTail && j < QSIZE; ++x){
      if (EMPTY(x)){
         ++i;
         continue;
      }
      
      extensionIndex[i] = extensionIndex[j];
      
      /* Sanity check: do the extension rows match the node rows? */
      if (extensionIndex[j] > 1){
         y = x;
         for (k = 0; k < 2*period; ++k){
            uint16_t startRow = extensionRows[extensionIndex[j]][1] + 1;
            if (extensionRows[extensionIndex[j]][startRow - k] != ROW(y)){
               fprintf(stderr, "Warning: non-matching rows detected at node %u in doCompactPart2()\n",x);
               free(extensionRows[extensionIndex[j]]);
               extensionRows[extensionIndex[j]] = 0;
               extensionIndex[i] = 0;
               break;
            }
            y = PARENT(y);
         }
      }
      if (j > i) extensionIndex[j] = 0;
      ++i;
      ++j;
   }
   for (j = qTail - qHead; j < QSIZE; ++j){
      extensionIndex[j] = 0;
   }
   extQHead = 0;
   extQTail = qTail - qHead;
}

void doCompact(void) {
   /* make sure we still have something left in the queue */
   if (qIsEmpty()) {
      qTail = qHead = 0;   /* nothing left, make an extremely compact queue */
      return;
   }
   /* First loop of part 1 requires qTail-1 to be non-empty.  Make it so. */
   while (EMPTY(qTail-1))
      qTail--;
   
   doCompactPart1();
   if (dumpFlag == DUMPPENDING) dumpState();
   doCompactPart2();
}

/* ================= */
/*  Lookahead cache  */
/* ================= */

#ifndef NOCACHE
int cacheLookup(uint64_t id1, uint64_t id2, uint64_t id3, int rowsKey) {
#ifndef QSIMPLE
   if (params[P_CACHEMEM] == 0) return 0;
#endif
   unsigned long long slot = (unsigned long long)id1 +
      17 * (unsigned long long)id2 + 257 * (unsigned long long)id3 +
      513ULL * (unsigned)rowsKey;
   slot ^= slot >> 16;
   slot = slot + (slot >> 15);
   slot &= (cachesize-1);
   cacheentry *entry = &(cache[omp_get_thread_num()][slot]);
   if (entry->id1 == id1 && entry->id2 == id2 && entry->id3 == id3 && entry->rowsKey == rowsKey)
      return -2 + entry->r;
   entry->id1 = id1;
   entry->id2 = id2;
   entry->id3 = id3;
   entry->rowsKey = rowsKey;
   return (int) slot;
}

void cacheStore(int slot, int result) {
#ifndef QSIMPLE
   if (params[P_CACHEMEM])
#endif
      cache[omp_get_thread_num()][slot].r = result;
}
#endif

/* ========================== */
/*  Primary search functions  */
/* ========================== */

void process(node theNode);
void bfsPrecompute(void);
extern int preN, preCursor;
int depthFirst(node theNode, uint16_t howDeep, const uint16_t **candEnd, int *candLeft, row *pathRows, _Atomic int *remainingItems, _Atomic int *forceExit, _Atomic int *passed);

static void deepen(void) {
   preN = preCursor = 0;   /* compaction moves nodes: drop precomputed lookaheads */
   /* compute amount to deepen, apply reduction if too deep */
   int deepeningAmount;
   int i = currentDepth();
   
   if (i >= lastDeep) deepeningAmount = MINDEEP;
   else deepeningAmount = lastDeep + MINDEEP - i;   /* go at least MINDEEP deeper */
   
   if (params[P_FIRSTDEEP]){
      deepeningAmount = params[P_FIRSTDEEP];
      params[P_FIRSTDEEP] = 0;
   }
   
   lastDeep = i + deepeningAmount;
   
   /* report what's happening */
   printf("%d, deepening %d, ", i, deepeningAmount);
   printCompactNumber(qTail - qHead);
   printf("/");
   printCompactNumber(qTail);
   fflush(stdout);
   
   _Atomic int remainingItems = 0;
   _Atomic int forceExit = 0;
   _Atomic int passed = 0;
   atomic_store_explicit(&remainingItems, qTail - qHead, memory_order_seq_cst);
   atomic_store_explicit(&forceExit, 0, memory_order_seq_cst);
   atomic_store_explicit(&passed, 0, memory_order_seq_cst);

   /* go through queue, deepening each one */
   #pragma omp parallel
   {
      const uint16_t **candEnd;
      int *candLeft;
      row *pathRows;
      
      candEnd = (const uint16_t**)calloc((deepeningAmount + 4 * params[P_PERIOD]), sizeof(*candEnd));
      candLeft = (int*)calloc((deepeningAmount + 4 * params[P_PERIOD]), sizeof(*candLeft));
      pathRows = (row*)calloc((deepeningAmount + 4 * params[P_PERIOD]), sizeof(*pathRows));
      if (!candEnd || !candLeft || !pathRows) {
         fprintf(stderr, "Error: out of memory while deepening\n");
         exit(1);
      }
      
      long long j;
      #pragma omp for schedule(dynamic, CHUNK_SIZE)
      for (j = qHead; j < qTail; j++) {
         if (!EMPTY(j) && !depthFirst((node)j, (uint16_t)deepeningAmount, candEnd, candLeft, pathRows, &remainingItems, &forceExit, &passed))
            MAKEEMPTY(j);
         atomic_fetch_sub_explicit(&remainingItems, 1, memory_order_relaxed);
      }
      free(candEnd);
      free(candLeft);
      free(pathRows);
   }
   
   /* before reporting new queue size, shrink tree back down */
   printf(" -> ");
   fflush(stdout);
   
   /* signal time for dump */
   if (params[P_DUMPMODE] != D_DISABLED && time(NULL) - lastDumpTime >= params[P_DUMPINTERVAL]){
      dumpFlag = DUMPPENDING;
      time(&lastDumpTime);
   }
   
   doCompact();
   
   /* now finish report */
   printCompactNumber(qTail - qHead);
   printf("/");
   printCompactNumber(qTail);
   printProgress();
   printf("\n");
   
   /* Report successful/unsuccessful dump */
   if (dumpFlag == DUMPSUCCESS) {
      timeStamp();
      printf("State dumped to %s\n",dumpFile);
      if (dumpNum == DUMPLIMIT){
         timeStamp();
         printf("Sequential dump limit reached.  Changing to overwrite mode.\n");
      }
   }
   else if (dumpFlag == DUMPFAILURE) {
      timeStamp();
      printf("State dump unsuccessful\n");
      if (dumpMode == D_OVERWRITE)
         dumpNum--;   /* reset dumpNum so that the next dump doesn't overwrite our "backup" */
   }
   dumpFlag = DUMPRESET;
   
   fflush(stdout);
}

static void breadthFirst(void) {
   while (!aborting && !qIsEmpty()){
      if (qTail - qHead >= (1LLU<<params[P_DEPTHLIMIT]) || qTail >= QSIZE - QSIZE/16){
         timeStamp();
         printf("Queue full, depth ");
         deepen();
      }
      else if (params[P_EVERYDEPTH] && qHead == nextRephase){
         timeStamp();
         printf("Depth ");
         deepen();
      }
      else
         { if (preCursor >= preN) bfsPrecompute(); process(dequeue()); }
   }
}

/* Save rows from deepening step in pathRows array (probably should be stored as a struct instead) */
void saveExtension(node theNode, uint32_t startRow, uint32_t howDeep, row *pathRows) {
   uint32_t extIndex;
   uint16_t totalDepth = (uint16_t) ((startRow + howDeep) & UINT16_MAX); /* truncate to fit in array */
   #pragma omp critical(findDeepIndex)
   {
      extIndex = 2;
      while (extIndex < 1LLU << (params[P_DEPTHLIMIT] + 1) && extensionRows[extIndex]) ++extIndex;
      if (extIndex == 1LLU << (params[P_DEPTHLIMIT] + 1)){
         fprintf(stderr,"Error: no available extension indices.\n");
         aborting = 1;
      }
      if (!aborting){
         extensionRows[extIndex] = (row*)calloc( totalDepth + 1 + 2,
                                                sizeof(**extensionRows) );
      }
   }
   if (aborting) return;
   
   memcpy( extensionRows[extIndex] + 2,
           pathRows,
           (totalDepth + 1) * (long long)sizeof(**extensionRows) );
   
   extensionRows[extIndex][0] = totalDepth;
   extensionRows[extIndex][1] = (uint16_t) startRow;
   
   extensionIndex[extQHead + theNode - qHead] = extIndex;
}

/* ========================== */
/*  Print usage instructions  */
/* ========================== */

void printHelp(void) {
   printf("Usage:    ./qfind "
#ifndef QSIMPLE
                             "-v <velocity> "
#endif
                                        "-w <width> -s <symmetry> [options...]\n"
          "       or\n"
          "          ./qfind -l <file> [options...]\n");
   printf("\n");
   printf("qfind is a program that searches for orthogonal spaceships and waves in Life\n"
          "and related cellular automata.  Options are read left to right, with subsequent\n"
          "occurrences of the same option overwriting the previous value.\n");
   printf("\n");
#ifdef QSIMPLE
   printf("You have compiled with the following values:\n");
   printf("\n");
   printf("Period: %d\n",PERIOD);
   printf("Offset: %d\n",OFFSET);
   printf("\n");
#endif
   printf("Required (except when loading from a saved state):\n");
#ifndef QSIMPLE
   printf("  -v, --velocity <velocity>     written in the form <translation>c/<period>\n");
#endif
   printf("  -w, --width <number>          logical width (full width depends on symmetry)\n");
   printf("  -s, --symmetry <(asymmetric|odd|even|gutter)>  spaceship symmetry type\n");
   printf("\n");
   printf("Search options:\n");
   printf("  -r, --rule <rule>             cellular automaton rule written in Hensel\n"
          "                                notation (Default: B3/S23)\n"
          "                                '~' is used to specify a list of forbidden\n"
          "                                conditions.  For example, -r B3~6c7/S23~8\n"
          "                                searches in B3/S23 for ships that never contain\n"
          "                                the B6c, B7, or S8 neighborhoods.\n");
#ifdef _OPENMP
   printf("  -t, --threads <number>        number of threads during deepening.  Default is\n"
          "                                one less than the number of available cores.\n");
#endif
   printf("  -f, --found <number>          maximum number of spaceships to output\n");
   printf("  -i, --increment <number>      minimum deepening increment (default: 3)\n");
   printf("  -g, --min-extension <number>  minimum length of saved extensions\n");
   printf("  -n, --first-depth <number>    depth of first deepening step\n");
   printf("      --fixed-depth <number>    deepen at every new depth by the given amount\n");
   printf("  -e, --extend <filename>       file containing the initial rows for a search.\n"
          "                                Use the Golly script get-rows.lua to easily\n"
          "                                generate the initial rows file.\n");
#ifdef _OPENMP
   printf("  (--enable-early-exit|--disable-early-exit)\n"
          "                                enable/disable early exit during deepening step\n"
          "                                when threads become idle (default: enabled)\n");
#endif
   printf("  --lookahead-depth <1|2|3>     how much of the lookahead test to apply\n"
          "                                (default: 3, the full test).  1 or 2 prune less\n"
          "                                but cost less per candidate.  The search stays\n"
          "                                exhaustive, though which variants of a ship get\n"
          "                                printed can differ (as with -i or -q).  Measured\n"
          "                                from 25%% faster to 60%% slower depending on the\n"
          "                                search.\n");
   printf("  --tables <auto|pair|factored> lookup table representation (default: auto =\n"
          "                                per-pair below width 10, factored from 10).\n"
          "                                factored needs a few MB instead of about 1.5 GB\n"
          "                                at width 10 (about 7x more per extra width) and\n"
          "                                is about as fast or faster from width 10; per-pair\n"
          "                                is faster below.  Exhaustive either way.\n");
   printf("  --lookahead-method <auto|cone|split|hybrid>\n"
          "                                how the full lookahead test is computed (default:\n"
          "                                auto = split from width 12, cone below).  All\n"
          "                                give the same answers, so the search is the same;\n"
          "                                only speed differs.  cone: whole rows (cost per\n"
          "                                test grows with width); split: low and high half\n"
          "                                rows joined on their 2 overlap bits (cost per test\n"
          "                                does not grow with width); hybrid: cone, switching\n"
          "                                to split when a test gets expensive.\n");
   printf("  --lookahead-budget <N>        hybrid: work units before switching to split\n"
          "                                (default: 256).\n");
   printf("  --split-cache-bits <N>        split and hybrid: 2^N cache entries per half,\n"
          "                                shared by all threads (default: %d = %d MB).\n"
          "                                Larger is faster at width 11 and above.\n",
          SPLIT_CACHE_BITS, (int)((288LL << SPLIT_CACHE_BITS) >> 20));
   printf("\n");
   printf("Memory options:\n");
   printf("  -c, --cache-mem <number>      allocate N megabytes per thread for a lookahead\n"
          "                                cache (default: 0 = disabled).  The cache uses the\n"
          "                                original lookahead code; the uncached fast\n"
          "                                lookahead measured 0-36%% faster at every speed\n"
          "                                tested, so caching is no longer on by default.\n");
   printf("  -m, --mem-limit <number>      limits lookup table memory to N megabytes\n");
   printf("  -q, --queue-bits <number>     set BFS queue size to 2^N nodes (default: %d)\n", QBITS);
   printf("  -h, --hash-bits <number>      set hash table size to 2^N nodes (default: %d)\n"
          "                                Use -h 0 to disable duplicate elimination.\n", HASHBITS);
   printf("  -b, --base-bits <number>      groups 2^N queue entries to an index node\n"
          "                                (default: 4)\n");
   printf("\n");
   printf("Save/load options:\n");
   printf("  -d, --dump-root <string>      dump filename prefix\n");
   printf("  -a, --dump-interval <number>  wait at least N seconds between dumps\n");
   printf("      --dump-mode <(overwrite|sequential|disabled)>\n"
          "                                set dump mode\n");
   printf("  -l, --load <filename>         load search state from the given dump file\n");
   printf("  -j, --split <number>          split loaded search state into at most N files\n");
   printf("  -p, --preview                 preview partial results from the loaded state\n");
   printf("\n");
   printf("Output options (enabled by default):\n");
#ifndef QSIMPLE
   printf("  (--enable-subperiod|--disable-subperiod)    enable/disable printing of\n"
          "                                              subperiodic results\n");
#endif
   printf("  (--enable-deep-print|--disable-deep-print)  enable/disable printing ships\n"
          "                                              during deepening step\n");
   printf("  (--enable-longest|--disable-longest)        enable/disable printing longest\n"
          "                                              partial result at end of search\n");
   printf("\n");
   printf("Wave options:\n");
   printf("  -o, --boundary-sym <(disabled|odd|even|gutter)>  boundary symmetry type for\n"
          "                                                   wave searches\n");
   printf("\n");
   printf("Documentation options:\n");
   printf("  --help                        print usage instructions and exit\n");
   printf("\n");
   printf("Performance notes (this build; details in CHANGES.md):\n");
   printf("  Progress lines end with [elapsed time, memory in use (peak), table memory].\n");
   printf("  Lookup tables use a compact layout (about 40%% less table memory).\n");
#ifdef PACKED
   printf("  This build was compiled with -DPACKLISTS: list entries are bit-packed to\n"
          "  the search width.  About 24%% less table memory again, about 25%% slower.\n");
#else
   printf("  For about 24%% less table memory (and about 25%% slower), recompile with\n"
          "  -DPACKLISTS.\n");
#endif
#ifdef __linux__
   printf("  Large tables are allocated in 2 MB-aligned blocks with transparent huge\n"
          "  pages requested; check /sys/kernel/mm/transparent_hugepage/enabled\n"
          "  (\"always\" or \"madvise\" lets this take effect).\n");
#endif
#ifndef QSIMPLE
   printf("  Lookahead caching is off by default at all speeds (it used to be on above\n"
          "  c/5): the per-level fast lookahead measured faster.  -c N turns it back on.\n");
#endif
   printf("  Compile with -DLA_VERIFY to cross-check the fast lookahead against the\n"
          "  original on every call (slow; aborts on any mismatch).\n");
#ifndef QSIMPLE
   printf("\n");
   printf("Example search:\n"
          "    ./qfind -v c/5 -w 9 -s even -r B3/S23 -t 2\n\n"
          "  Searches Life (rule B3/S23) for c/5 orthogonal spaceships with even\n"
          "  bilateral symmetry and logical width 9 (full width 18) using two threads.\n\n");
#endif
   exit(0);
}

/* ======================== */
/*  Echo loaded parameters  */
/* ======================== */

void echoParams(void) {
   printf("\n");
   printf("Rule: %s\n",rule);
   printf("speed: ");
   if (params[P_OFFSET] != 1) printf("%d",params[P_OFFSET]);
   printf("c/%d\n", params[P_PERIOD]);
   printf("Width: %d\n", params[P_WIDTH]);
   if (params[P_SYMMETRY] == SYM_ASYM) printf("Symmetry: asymmetric\n");
   else if (params[P_SYMMETRY] == SYM_ODD) printf("Symmetry: odd\n");
   else if (params[P_SYMMETRY] == SYM_EVEN) printf("Symmetry: even\n");
   else if (params[P_SYMMETRY] == SYM_GUTTER) printf("Symmetry: gutter\n");
   if (params[P_BOUNDARYSYM] != SYM_UNDEF){
      printf("Wave search enabled\nBoundary symmetry: ");
      if (params[P_BOUNDARYSYM] == SYM_ODD) printf("odd\n");
      else if (params[P_BOUNDARYSYM] == SYM_EVEN) printf("even\n");
      else if (params[P_BOUNDARYSYM] == SYM_GUTTER) printf("gutter\n");
   }
#ifndef QSIMPLE
   if (params[P_FULLPERIOD] && gcd(period,offset)>1) printf("Suppress subperiodic results\n");
#endif
   if (params[P_DUMPMODE] != D_DISABLED){
      printf("Dump interval: %d second%s\n", params[P_DUMPINTERVAL], params[P_DUMPINTERVAL] == 1 ? "" : "s");
      printf("Dump mode: %s\n", params[P_DUMPMODE] == D_OVERWRITE ? "overwrite" : "sequential");
   }
   else
      printf("Dump disabled\n");
   printf("Queue size: 2^%d\n",params[P_QBITS]);
   printf("Hash table size: 2^%d\n",params[P_HASHBITS]);
   if (params[P_EVERYDEPTH])
      printf("Fixed deepening amount: %ld\n",
               params[P_FIRSTDEEP] ? (long)params[P_FIRSTDEEP] : lastDeep - currentDepth());
   else
      printf("Minimum deepening increment: %d\n",MINDEEP);
   if (params[P_PRINTDEEP] == 0) printf("Output disabled while deepening\n");
#ifndef NOCACHE
   if (params[P_CACHEMEM])
      printf("Cache memory per thread: %d megabytes\n", params[P_CACHEMEM]);
   else
      printf("Lookahead caching disabled\n");
#endif
   if (params[P_MEMLIMIT] >= 0) printf("Memory limit: %d megabytes\n",params[P_MEMLIMIT]);
#ifdef _OPENMP
   printf("Number of threads: %d\n",params[P_NUMTHREADS]);
#endif
   if (params[P_MINEXTENSION]) printf("Save depth-first extensions of length at least %d\n",params[P_MINEXTENSION]);
   if (params[P_LONGEST] == 0) printf("Printing of longest partial result disabled\n");
   if (lookaheadDepth != 3) printf("Lookahead depth: %d\n", lookaheadDepth);
   printf("Lookup tables: %s\n", useFactored ? "factored (two global half-tables)" : "one per row pair");
   if (lookaheadDepth == 3) printf("Lookahead method: %s\n", useSplitLookahead == 1 ? "split (half-row cones)" : useSplitLookahead == 2 ? "hybrid (cone, split when expensive)" : "cone (full rows)");
   if (lookaheadDepth == 3 && useSplitLookahead) printf("Split lookahead cache: %lld MB, shared by all threads (--split-cache-bits %d)\n", (288LL << splitCacheBits) >> 20, splitCacheBits);
   printf("\n");
}

/* ========================= */
/*  Preview partial results  */
/* ========================= */

static void preview(void /*int allPhases*/) {
   node i,j;
   row *pathRows;
   // int ph;  /* used for allPhases option */
   // node k;  /* used for allPhases option */

   for (i = qHead; (i<qTail) && EMPTY(i); i++);
   for (j = qTail-1; (j>i) && EMPTY(j); j--);
   if (j<i) return;
   
   while (j>=i && !aborting) {
      if (!EMPTY(j)) {
         uint32_t extIndex = extensionIndex[extQHead + j - qHead];
         
         if (extIndex > 1){
            pathRows = (row*) malloc((2*period + 1 + extensionRows[extIndex][0]
                                   - extensionRows[extIndex][1] + 1) * (long long)sizeof(*pathRows));
            int m;
            node x = j;
            for (m = 2*period; m >= 0; --m){
               pathRows[m] = ROW(x);
               x = PARENT(x);
            }
            memcpy(pathRows + 2*period+1,
                   extensionRows[extIndex] + 2 + extensionRows[extIndex][1],
                   (extensionRows[extIndex][0] - extensionRows[extIndex][1] + 1) * (long long)sizeof(*pathRows));
            uint32_t currRow = 2*period + 1  + extensionRows[extIndex][0]
                                             - extensionRows[extIndex][1]
                                             + 1;
            success(j, pathRows, 2*period, currRow - 1);
            free(pathRows);
         }
         else{
            success(j, NULL, 0, 0);
         }
         /*
         if (allPhases) {
            k=j;
            for (ph = 1; ph < period; ph++) {
               k=PARENT(k);
               success(k, NULL, 0, 0);
            }
         }
         */
      }
      j--;
   }
}

/* =============================== */
/*  Check parameters for validity  */
/* =============================== */

int loadDumpFlag = 0;
int previewFlag = 0;
int initRowsFlag = 0;

void optError(const char *errorMsg, const char *opt) {
   fprintf(stderr, "Error: %s%s\n", errorMsg, opt);
   aborting = 1;
}

void printError(const char *errorMsg) {
   optError(errorMsg,"");
}

/* Reads and reports values from ruleTable for a range of conditions.
** input:  a string of birth conditions or a string
**         of survival conditions (not both)
**         (e.g., "B34-w6ci" or "S0123")
**                                                   
** output: the function returns -1 if all input
**         conditions have value -1 in ruleTable;
**         otherwise, it returns the value in
**         ruleTable of the input conditions if they
**         are all the same (ignoring -1), and it
**         returns 2 if they are not the same.
*/
int checkConditions(const char *p) {
   int letterHits[256];
   int ruleOffset = 0;
   int i;
   if (*p == 's' || *p == 'S') ruleOffset = 256;
   p++;
   int val = -1;
   while (*p != '\0'){
      char digit = *p++;
      int negCount = 0;
      if (*p == '\0' || ('0' <= *p && *p <= '8')){
         for (i = 0; i < 256; i++){
            if (rulekeys[i][0] == digit && ruleTable[ruleOffset+i] != -1){
               if (val == -1)
                  val = ruleTable[ruleOffset+i];
               if (val != ruleTable[ruleOffset+i])
                  return 2;
            }
         }
      }
      else if (*p == '-'){
         p++;
         for (i = 0; i < 256; i++)
            letterHits[i] = 0;
         for (; *p != '\0' && !('0' <= *p && *p <= '8'); p++){
            for (i = 1; i < 256; i++)
               if (rulekeys[i][0] == digit && *p != rulekeys[i][1])
                  letterHits[i]++;
            negCount++;
         }
         for (i = 0; i < 256; i++){
            if (letterHits[i] == negCount && ruleTable[ruleOffset+i] != -1){
               if (val == -1)
                  val = ruleTable[ruleOffset+i];
               if (val != ruleTable[ruleOffset+i])
                  return 2;
            }
         }
      }
      else {
         for (; *p != '\0' && !('0' <= *p && *p <= '8'); p++){
            for (i = 0; i < 256; i++){
               if (rulekeys[i][0] == digit && rulekeys[i][1] == *p && ruleTable[ruleOffset+i] != -1){
                  if (val == -1)
                     val = ruleTable[ruleOffset+i];
                  if (val != ruleTable[ruleOffset+i])
                     return 2;
               }
            }
         }
      }
   }
   return val;
}

void checkRule(void) {
   /* Errors: no meaningful results possible */
   if (checkConditions("B0") == 1)
      printError("rules with B0 are not supported.");
   
   /* Included below are conditions that prevent spaceships from existing:
   ** For proofs of these conditions, see here:
   ** https://conwaylife.com/forums/viewtopic.php?f=11&t=5471
   ** https://ics.uci.edu/~eppstein/ca/glider.c
   **
   **       x <= 0: none of the conditions are satisfied in
   **               the minimum rule
   ** (x+1)%2 == 0: all of the conditions are satisfied in
   **               the maximum rule
   */
   if (checkConditions("B0") == -1){
      printError("any pattern that is not infinite in both dimensions must contain the B0\n       "
                 "neighborhood.");
   }
   if (checkConditions("B1c") == -1){
      printError("spaceships and waves must contain the B1c neighborhood.");
   }
   else if (checkConditions("B1e2a") == -1) {
      printError("spaceships and waves must contain at least one of the B1e or B2a\n       "
                 "neighborhoods.");
   }
   if (checkConditions("B1c") == 1 && checkConditions("B0") == 0){
      printError("patterns in rules with B1c and without B0 expand in all directions.");
   }
   else if (checkConditions("B1e2a") == 1 && checkConditions("B0") == 0) {
      printError("patterns in rules with B1e2a and without B0 expand in all directions.");
   }
   /* The following checks are done only for spaceship searches */
   if ( (params[P_BOUNDARYSYM] == SYM_UNDEF || params[P_SYMMETRY] == SYM_ASYM) ){ 
      if (checkConditions("B012ac3i") <= 0){
         printError("patterns in rules without any of B012ac3i cannot leave their initial\n       "
                    "bounding box.");
      }
      if (checkConditions("B012ae3a") <= 0) {
         printError("patterns in rules without any of B012ae3a cannot leave their initial\n       "
                    "bounding diamond.");
      }
      if (checkConditions("B01245") <=0 && checkConditions("S012345") <= 0) {
         printError("patterns in rules without any of B01245/S012345 cannot move a distance\n       "
                    "of more than one cell outside their initial bounding diamond.");
      }
      if ( checkConditions("B01e2a") <= 0 && 2 * params[P_OFFSET] > params[P_PERIOD]
                                          && params[P_PERIOD] > 0 ){
         printError("orthogonal spaceship speed limit in rules without any of B01e2a is c/2.");
      }
   }
   
   /* Warnings: no spaceships exist, but maybe we can get some interesting wickstretchers */
   if ( (params[P_BOUNDARYSYM] == SYM_UNDEF || params[P_SYMMETRY] == SYM_ASYM) ){
      if (checkConditions("B0") == 0 && (checkConditions("B23")+1)%2 == 0 && (checkConditions("S0")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of B23/S0 and without B0,\n"
                         "         because the trailing edge of a pattern cannot die.\n");
      }
      else if (checkConditions("B0") == 0 && checkConditions("B123") >= 1 && (checkConditions("S0123")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with one of B1, B2, or B3, all of S0123,\n"
                         "         and without B0, because the trailing edge of a pattern cannot die.\n");
      }
      if ((checkConditions("S012acek3aijn4a")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of S012acek3aijn4a and\n"
                         "         without B0, because patterns cannot shrink.\n");
      }
      if ((checkConditions("S1234-wz5-aqr6ce")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of S1234-wz5-aqr6ce and\n"
                         "         without B0, because connected patterns cannot shrink.\n");
      }
      if ((checkConditions("B34")+1)%2 == 0 && (checkConditions("S12345")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of B34/S12345 and without B0,\n"
                         "         because connected patterns cannot shrink.\n");
      }
      if ((checkConditions("B345")+1)%2 == 0 && (checkConditions("S1234")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of B345/S1234 and without B0,\n"
                         "         because connected patterns cannot shrink.\n");
      }
      if (checkConditions("B012") <= 0 && (checkConditions("S234567")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of S234567 and none of B012,\n"
                         "         because patterns cannot escape their bounding diamond without an\n"
                         "         immortal triangle.\n");
      }
   }
}

void checkGutter(void) {
   int i = 0;
   /* if there are forbidden birth conditions, i will be less than 256 */
   while (i < 256 && ruleTable[i] != -1)
      i++;
   
   if (checkConditions("B2ci4ci6i") <= 0)
      gutterSkew = 0;
   else if (checkConditions("B1c2kn3ny4yz5r6i") <= 0)
      gutterSkew = 1;
   else if (checkConditions("B12aikn3cqr4cnyz5er6i") <= 0)
      gutterSkew = 2;
   else {
      printError("gutters do not work with the given birth conditions.\n       "
                 "The forbidden birth conditions for different gutter types are\n       "
                 "  Skew 0: B2ci4ci6i\n       "
                 "  Skew 1: B1c2kn3ny4yz5r6i\n       "
                 "  Skew 2: B12aikn3cqr4cnyz5er6i");
   }
   if (gutterSkew && i < 256)
      fprintf(stderr, "Warning: forbidden birth conditions cannot be checked along a skew gutter.\n");
}

void checkParams(void) {
   const char *ruleError;
   
   /* Errors */
   
   ruleError = parseRule(rule, ruleTable);
   
   if (ruleError != 0){
      optError("failed to parse rule ", rule);
      fprintf(stderr, "       %s\n", ruleError);
   }
   else
      checkRule();
   
   if (params[P_SYMMETRY] == SYM_GUTTER || params[P_BOUNDARYSYM] == SYM_GUTTER)
      checkGutter();
   
#ifdef QSIMPLE
   if (gcd(PERIOD,OFFSET) > 1)
      printError("gcd(PERIOD,OFFSET) > 1 not supported when compiled with preset values.\n       "
                 "Please recompile with PERIOD and OFFSET undefined.");
#else
   if (params[P_PERIOD] > MAXPERIOD)
      printError("maximum allowed period (" XSTR(MAXPERIOD) ") exceeded.");
   if (params[P_OFFSET] > params[P_PERIOD] && params[P_PERIOD] > 0)
      printError("translation cannot exceed period.");
   if (params[P_OFFSET] == params[P_PERIOD] && params[P_PERIOD] > 0)
      printError("photon searches are not supported.");
#endif
   if (params[P_PERIOD] == 0)
      printError("you must specify a velocity (-v).");
   if (params[P_WIDTH] == 0)
      printError("you must specify a width (-w).");
   if (params[P_SYMMETRY] == SYM_UNDEF)
      printError("you must specify a symmetry type (-s).");
   if (params[P_BOUNDARYSYM] == SYM_ASYM)
      printError("asymmetric wave searching is not supported.");
   if (previewFlag && !loadDumpFlag)
      printError("the search state must be loaded from a file to preview partial results.\n");
   if (initRowsFlag && loadDumpFlag){
      printError("initial rows file cannot be used when the search state is loaded from a\n       "
                 "saved state.");
   }
   if (params[P_WIDTH] > MAXWIDTH)
      printError("maximum allowed width (14) exceeded.");
   if (params[P_NUMTHREADS] <= 0)
      printError("number of threads (-t) must be positive.");
   if (params[P_QBITS] < 4)
      printError("queue bits (-q) must be at least 4.");
   if (params[P_MINDEEP] > 10000 || params[P_FIRSTDEEP] > 10000)
      printError("deepening amounts (-i, -n, --fixed-depth) must not exceed 10000.");
   if (params[P_BASEBITS] <= 0)
      printError("base bits (-b) must be positive.");
   if (params[P_BASEBITS] >= params[P_QBITS])
      printError("base bits (-b) must be less than queue bits (-q).");
   if (params[P_HASHBITS] < 0)
      printError("hash bits (-h) must be nonnegative.");
   
   /* Warnings */
   if (2 * params[P_OFFSET] > params[P_PERIOD] && params[P_PERIOD] > 0){
      fprintf(stderr, "Warning: searches for speeds exceeding c/2 may not work correctly.\n");
   }
#ifndef NOCACHE
   if (params[P_CACHEMEM] > 0){
      fprintf(stderr, "Note: lookahead caching (-c) uses the original lookahead, which measured\n"
                      "      slower than the uncached fast lookahead.  Use -c 0 to disable it.\n");
      if (lookaheadDepth != 3)
         fprintf(stderr, "Warning: --lookahead-depth is ignored while lookahead caching is enabled.\n");
   }
#endif
   if (params[P_SYMMETRY] == SYM_ASYM && params[P_BOUNDARYSYM] != SYM_UNDEF){
      fprintf(stderr, "Warning: the wave symmetry settings are equivalent to a spaceship search.\n");
      params[P_SYMMETRY] = params[P_BOUNDARYSYM];
      params[P_BOUNDARYSYM] = SYM_UNDEF;
   }
   /* Reduce values to prevent integer overflow */
   if (params[P_QBITS] > 31 && !aborting){
      fprintf(stderr, "Warning: queue bits (-q) reduced to 31.\n");
      params[P_QBITS] = 31;      /* corresponds to a queue size of 2GB */
      if (params[P_BASEBITS] >= params[P_QBITS]){
         fprintf(stderr, "Warning: base bits (-b) reduced to 30.\n");
         params[P_BASEBITS] = 30;
      }
   }
   if (params[P_HASHBITS] > 31 && !aborting){
      fprintf(stderr, "Warning: hash bits (-h) reduced to 31.\n");
      params[P_HASHBITS] = 31;   /* corresponds to a hash table size of 2GB */
   }
   
}

/* ============================ */
/*  Load saved state from file  */
/* ============================ */

char * loadFile;

void loadFail(void) {
   fprintf(stderr, "Load from file %s failed\n", loadFile);
   exit(1);
}

signed int loadInt(FILE *fp) {
   signed int v;
   if (fscanf(fp,"%d\n",&v) != 1) loadFail();
   return v;
}

unsigned long loadUInt(FILE *fp) {
   unsigned long v;
   if (fscanf(fp,"%lu\n",&v) != 1) loadFail();
   return v;
}

void loadParams(void) {
   FILE * fp;
   unsigned int i;
   
   fp = fopen(loadFile, "r");
   if (!fp) loadFail();
   if (loadUInt(fp) != FILEVERSION){
      printf("Incompatible file version\n");
      exit(1);
   }
   
   /* Load rule */
   if (fscanf(fp,"%150s\n",loadRule) != 1) loadFail();
   rule = loadRule;
   
   /* Load dump root */
   if (fscanf(fp,"%250s\n",loadDumpRoot) != 1) loadFail();
   dumpRoot = loadDumpRoot;
   
   /* Load parameters */
   for (i = 0; i < NUM_PARAMS; ++i)
      params[i] = loadInt(fp);
}

void loadState(void) {
   FILE * fp;
   unsigned long i, j;
   int k;
   
   fp = fopen(loadFile, "r");
   if (!fp) loadFail();
   
   /* Skip lines that are loaded in loadParams() */
   loadUInt(fp);                                      /* skip file version */
   if (fscanf(fp, "%*[^\n]\n") != 0) loadFail();      /* skip rule */
   if (fscanf(fp, "%*[^\n]\n") != 0) loadFail();      /* skip dump root */
   for (i = 0; i < NUM_PARAMS; ++i) loadInt(fp);      /* skip parameters */
   
   /* Load / initialise globals */
   width          = loadInt(fp);
   period         = loadInt(fp);
   offset         = loadInt(fp);
   lastDeep       = loadInt(fp);
   if (!splitNum)
      dumpNum     = loadInt(fp);
   else {
      if (fscanf(fp, "%*[^\n]\n") != 0) loadFail();   /* skip dumpNum when splitting loaded state */
   }
   if (params[P_DUMPMODE] == D_SEQUENTIAL)
      dumpNum = 1;
   
   
   aborting        = 0;
   nRowsInState    = period+period;   /* how many rows needed to compute successor graph? */
   
   params[P_DEPTHLIMIT] = DEFAULT_DEPTHLIMIT;
   
    /* Allocate space for the data structures */
   queueBase = (node*)malloc((QSIZE>>BASEBITS)*sizeof(node));
   queueRows = (row*)malloc(QSIZE*sizeof(row));
   if (queueBase == 0 || queueRows == 0) {
      printf("Unable to allocate BFS queue!\n");
      exit(1);
   }
   
   if (hashBits == 0) seenHash = 0;
   else {
      seenHash = (node*)malloc(HASHSIZE*sizeof(node));
      if (seenHash == 0) printf("Unable to allocate hash table, duplicate elimination disabled\n");
   }
   
   /* Load up BFS queue */
   qHead  = (node) loadUInt(fp);
   qEnd   = (node) loadUInt(fp);
   qStart = QSIZE - qEnd;
   qEnd   = QSIZE;
   qHead += qStart;
   if (qStart > QSIZE || qStart < QSIZE/16 || qHead > qEnd) {
      printf("BFS queue is too small for saved state\n");
      exit(1);
   }
   for (i = qStart; i < qEnd; ++i)
      queueRows[i] = (row) loadUInt(fp);
   
   /* Load extension rows for each queue node */
   extensionRows = (row**)calloc(1LLU << (params[P_DEPTHLIMIT] + 1), sizeof(*extensionRows));
   extensionIndex = (uint32_t*)calloc(EXTINDEX_SIZE, sizeof(*extensionIndex));
   
   uint32_t extIndex = 2;
   extQTail = 0;
   
   while ( (k = fscanf(fp,"%lu\n",&j)) != EOF){
      if (k == 0)
         loadFail();
      if (j == 0){
         j = loadUInt(fp);
         if (j > QSIZE - extQTail) loadFail();
         for (i = 0; i < j; ++i){
            extensionIndex[extQTail] = 1;
            ++extQTail;
         }
         continue;
      }
      if (extQTail >= QSIZE || extIndex >= 1LLU << (params[P_DEPTHLIMIT] + 1) || j > UINT16_MAX) loadFail();
      extensionRows[extIndex] = (row*)calloc( j + 1 + 2, sizeof(**extensionRows));
      extensionRows[extIndex][0] = (row) j;
      for (i = 1; i < j + 1 + 2; ++i){
         extensionRows[extIndex][i] = (row) loadUInt(fp);
      }
      extensionIndex[extQTail] = extIndex;
      ++extIndex;
      ++extQTail;
   }
   
   fclose(fp);

   flipRow = (row*)malloc(sizeof(*flipRow)<<width);
   makeFlip();
   
   /* complete compaction */
   doCompactPart2();
   
   /* Let the user know that we got this far (suppress if splitting) */
   if (!splitNum) printf("State successfully loaded from file %s\n",loadFile);
   
   fflush(stdout);
}

/* ================================================= */
/*  Load initial rows for extending partial results  */
/* ================================================= */

void printRow(row theRow) {
   int i;
   for (i = width - 1; i >= 0; --i) printf("%c",(theRow & 1 << i ? 'o' : '.'));
   printf("\n");
}

void loadInitRows(char * file) {
   FILE * fp;
   int i,j;
   char rowStr[MAXWIDTH + 2];
   row theRow = 0;
   
   loadFile = file;
   fp = fopen(loadFile, "r");
   if (!fp) loadFail();
   
   printf("Starting search from rows in %s:\n",loadFile);
   
   for (i = 0; i < 2 * period; i++){
      if (fscanf(fp,"%15s",rowStr) != 1 || (int)strlen(rowStr) < width) loadFail();
      for (j = 0; j < width; j++){
         theRow |= ((rowStr[width - j - 1] == '.') ? 0:1) << j;
      }
      printRow(theRow);
      enqueue(dequeue(),theRow);
      theRow = 0;
   }
   fclose(fp);
}

/* ============== */
/*  Set Defaults  */
/* ============== */

void setDefaultParams(void) {
#ifdef QSIMPLE
   params[P_PERIOD] = PERIOD;
   params[P_OFFSET] = OFFSET;
#else
   params[P_PERIOD] = 0;
   params[P_OFFSET] = 0;
#endif
   params[P_WIDTH] = 0;
   params[P_SYMMETRY] = SYM_UNDEF;
   params[P_REORDER] = 1;     /* 0 and 2 are also valid values, but this cannot currently be set at runtime. */
   params[P_DUMPINTERVAL] = 1800;    /* 30 minutes */
   params[P_BASEBITS] = 4;
   params[P_QBITS] = QBITS;
   params[P_HASHBITS] = HASHBITS;
   params[P_NUMTHREADS] = omp_get_num_procs() > 1 ? omp_get_num_procs() - 1 : 1;
   params[P_MINDEEP] = 3;
   /* A negative value for params[P_CACHEMEM] means use that amount of  */
   /* memory if speed > c/5 and turn off caching otherwise.  A positive */
   /* value forces caching even if speed <= c/5.                        */
   params[P_CACHEMEM] = 0;   /* off: the fast lookahead beats the cache at every speed tested */
   params[P_MEMLIMIT] = -1;   /* negative value means no limit */
   params[P_PRINTDEEP] = 1;
   params[P_LONGEST] = 1;
   params[P_FIRSTDEEP] = 0;
   params[P_NUMSHIPS] = 0;
   params[P_MINEXTENSION] = 0;
   params[P_FULLPERIOD] = 0;
   params[P_BOUNDARYSYM] = SYM_UNDEF;
   params[P_DUMPMODE] = D_OVERWRITE;
   params[P_EVERYDEPTH] = 0;
   params[P_EARLYEXIT] = 1;
}

/* =============== */
/*  Parse options  */
/* =============== */

#define no_argument       (0)
#define required_argument (1)
#define optional_argument (2)

struct option {
   const char *name;
   int has_arg;
   int val;
};

/* Nonstandard getopt_long():
** This is a simple implementation of a getopt-like function.
** It's missing some features that ordinary getopt_long() has,
** but it has mostly the same input and output.
**
** Differences from getopt_long():
**   I/O:   char **optName and char **optArg:
**             used to save pointers to option names and arguments.
**             *optName and *optArg will point to elements of argv[].
**             If no argument is provided, *optArg will be 0.
**             If my_getopt() returns -1, then *optName and *optArg
**             will be unchanged.
*/
int my_getopt( int argc,
               char *argv[],
               const char *shortOpts,
               struct option *longOpts,
               char **optName,
               char **optArg   )
{
   static int i = 0;
   const char *charInd;
   if (++i == argc)  /* we're at the end of argv */
      return -1;
   *optName = argv[i];
   *optArg  = NULL;
   if (argv[i][0] == '-'){
      if (argv[i][1] == '\0')
         return '?';
      /* check if short option is contained in short options list */
      else if (argv[i][1] != '-' && argv[i][2] == '\0' && (charInd = strchr(shortOpts, argv[i][1]))){
         if (*(charInd+1) == ':'){
            if (argc == i+1){  /* no argument */
               if (*(charInd+2) != ':')
                  return (shortOpts[0] == ':' ? ':' : '?');
            }
            else{
               *optArg = argv[++i];
               return (int)argv[i-1][1];
            }
         }
         return (int)argv[i][1];
      }
      /* check if long option is contained in long options list */
      else if (argv[i][1] == '-') {
         while (longOpts->name != 0 && strcmp(argv[i] + 2, longOpts->name)) 
            longOpts++;
         if (longOpts->name){
            if (longOpts->has_arg){
               if (argc == i+1){  /* no argument */
                  if (longOpts->has_arg != optional_argument)
                     return (shortOpts[0] == ':' ? ':' : '?');
               }
               else
                  *optArg = argv[++i];
            }
            return longOpts->val;
         }
      }
   }
   return '?';
}

int readInt(char *opt, char *arg) {
   int i = 0;
   char c;
   if (arg == 0)
      aborting = 1;
   else if (strlen(arg) > 9){    /* prevent integer overflow in sscanf() from bad input */
      fprintf(stderr, "Error: invalid argument %s in option %s.\n", arg, opt);
      aborting = 1;
   }
   else if (sscanf(arg, "%d%c", &i, &c) != 1){
      fprintf(stderr, "Error: invalid argument %s in option %s.\n", arg, opt);
      aborting = 1;
   }
   return i;
}

const char *parseVelocity(char *velString, int *per, int *off) {
   int xoff = 0;
   *per = 1;
   *off = 1;
   char b = 0;
   char c = 0;
   if (strlen(velString) > 11)   /* prevent integer overflow in sscanf() from bad input */
      return "velocity string too long";
   else if (!strcmp(velString,"c"))
      return 0;
   else if (sscanf(velString, "c/%d%c%c", per, &b, &c) >= 1){
      if (b == 'd' && !c)
         return "diagonal spaceship searches are not supported.";
      else if (b == '\0' || (b == 'o' && !c))
         return 0;
      else
         return "illegal characters after velocity";
   }
   else if (sscanf(velString, "%dc/%d%c%c", off, per, &b, &c) >= 2){
      if (*off == 0)
         return "oscillator searches are not supported.";
      else if (b == 'd' && !c)
         return "diagonal spaceship searches are not supported.";
      else if (*off < 0)
         return "offset must be positive.";
      else if (b == '\0' || (b == 'o' && !c))
         return 0;
      else
         return "illegal characters after velocity";
   }
   else if (sscanf(velString, "(%d,%d)c/%d%c", off, &xoff, per, &c) >= 3){
      if (c != '\0' && c != 'o')
         return "illegal characters after velocity";
      else if (xoff != 0) {
         if (*off == 0){
            *off = xoff;
            xoff = 0;
            return (*off < 0) ? "offset must be positive." : 0;
         }
         else if (xoff == *off || xoff * (-1) == *off)
            return "diagonal spaceship searches are not supported.";
         else
            return "oblique spaceship searches are not supported.";
      }
      else if (*off == 0)
         return "oscillator searches are not supported.";
      else if (*off < 0)
         return "offset must be positive.";
      else
         return 0;
   }
   return "Unable to read offset and period.";
}

void parseOptions(int argc, char *argv[]) {
   char *optArg = 0;
   char *optName = 0;
   int c;
   
   if (argc <= 1){
      printf("\n");
      printHelp();
   }
   printf("Input:");
   for (c = 1; c < argc; c++)
      printf(" %s", argv[c]);
   printf("\n\n");
   
   /* list of long options */
   struct option options[] = {
   // {"",                    no_argument,        -1},   /* end option parsing when "--" encountered */
      {"help",                no_argument,       256},
      {"rule",                required_argument, 'r'},
      {"width",               required_argument, 'w'},
      {"symmetry",            required_argument, 's'},
      {"boundary-sym",        required_argument, 'o'},
      {"boundary-symmetry",   required_argument, 'o'},
      {"mem-limit",           required_argument, 'm'},
      {"memory-limit",        required_argument, 'm'},
      {"cache-mem",           required_argument, 'c'},
      {"cache-memory",        required_argument, 'c'},
      {"first-depth",         required_argument, 'n'},
      {"increment",           required_argument, 'i'},
      {"queue-bits",          required_argument, 'q'},
      {"hash-bits",           required_argument, 'h'},
      {"base-bits",           required_argument, 'b'},
#ifdef _OPENMP
      {"threads",             required_argument, 't'},
#endif
      {"found",               required_argument, 'f'},
      {"min-extension",       required_argument, 'g'},
      {"minimum-extension",   required_argument, 'g'},
      {"extend",              required_argument, 'e'},
      {"dump-root",           required_argument, 'd'},
      {"load",                required_argument, 'l'},
      {"split",               required_argument, 'j'},
      {"dump-interval",       required_argument, 'a'},
      {"dump-int",            required_argument, 'a'},
      {"preview",             no_argument,       'p'},
#ifndef QSIMPLE
      {"velocity",            required_argument, 'v'},
      {"enable-subperiod",    no_argument,       257},
      {"enable-subperiodic",  no_argument,       257},
      {"disable-subperiod",   no_argument,       258},
      {"disable-subperiodic", no_argument,       258},
#endif
      {"enable-deep-print",   no_argument,       259},
      {"disable-deep-print",  no_argument,       260},
      {"enable-longest",      no_argument,       261},
      {"disable-longest",     no_argument,       262},
      {"dump-mode",           required_argument, 263},
      {"fixed-depth",         required_argument, 264},
#ifdef _OPENMP
      {"enable-early-exit",   no_argument,       265},
      {"disable-early-exit",  no_argument,       266},
#endif
      {"lookahead-depth",     required_argument, 267},
      {"tables",              required_argument, 268},
      {"lookahead-method",    required_argument, 269},
      {"lookahead-budget",    required_argument, 270},
      {"split-cache-bits",    required_argument, 271},
      {0, 0, 0}   /* marks end of long options list */
   };
   
   while ( (c = my_getopt( argc,
                           argv,
                           ":"   /* List of short options; "<option>:" means argument required */
#ifndef QSIMPLE
                           "v:V:"
#endif
#ifdef _OPENMP
                           "t:T:"
#endif
                           "a:b:c:d:e:f:g:h:i:j:l:m:n:o:pq:r:s:w:"     /* Currently unused:   */
                           "A:B:C:D:E:F:G:H:I:J:L:M:N:O:PQ:R:S:W:",    /* k,u,x,y,z,K,U,X,Y,Z */
                           options,
                           &optName,
                           &optArg)) != -1 )
   {
      switch (c) {
         case 'r': case 'R':
            rule = optArg;
            if (strlen(rule) > 150)  /* any valid rule can be written in 141 characters or fewer */
               printError("rule string exceeds maximum allowed length (150).\n       "
                          "You must write the rule more efficiently.\n");
            break;
#ifndef QSIMPLE
         case 'v': case 'V':
         {
            const char *velError = 0;
            velError = parseVelocity(optArg, &params[P_PERIOD], &params[P_OFFSET]);
            if (velError){
               optError("invalid velocity ", optArg);
               fprintf(stderr, "       %s\n", velError);
               /* prevent additional error message in checkParams() */
               params[P_PERIOD] = 2; params[P_OFFSET] = 1;
            }
            if (params[P_PERIOD] <= 0){
               optError("invalid velocity ", optArg);
               fprintf(stderr, "       Period must be positive\n");
               /* prevent additional error message in checkParams() */
               params[P_PERIOD] = 2; params[P_OFFSET] = 1;
            }
            break;
         }
         case 257:   /* --enable-subperiod */
            params[P_FULLPERIOD] = 0;
            break;
         case 258:   /* --disable-subperiod */
            params[P_FULLPERIOD] = 1;
            break;
#endif
         case 'w': case 'W':
            params[P_WIDTH] = readInt(optName, optArg);
            if (params[P_WIDTH] <= 0){
               printError("width must be positive");
               params[P_WIDTH] = 1;    /* prevent additional error message in checkParams() */
            }
            break;
         case 's': case 'S':
            switch(optArg[0]) {
               case 'a': case 'A':
                  params[P_SYMMETRY] = SYM_ASYM; break;
               case 'o': case 'O':
                  params[P_SYMMETRY] = SYM_ODD; break;
               case 'e': case 'E':
                  params[P_SYMMETRY] = SYM_EVEN; break;
               case 'g': case 'G':
                  params[P_SYMMETRY] = SYM_GUTTER; break;
               default:
                  optError("unrecognized symmetry type ", optArg);
                  break;
            }
            break;
         case 'o': case 'O':
            switch(optArg[0]) {
               case 'a': case 'A':
                  params[P_BOUNDARYSYM] = SYM_ASYM; break;
               case 'o': case 'O':
                  params[P_BOUNDARYSYM] = SYM_ODD; break;
               case 'e': case 'E':
                  params[P_BOUNDARYSYM] = SYM_EVEN; break;
               case 'g': case 'G':
                  params[P_BOUNDARYSYM] = SYM_GUTTER; break;
               case 'd': case 'D':
                  params[P_BOUNDARYSYM] = SYM_UNDEF; break;
               default:
                  optError("unrecognized symmetry type ", optArg);
                  break;
            }
            break;
         case 'm': case 'M':
            params[P_MEMLIMIT] = readInt(optName, optArg);
            break;
         case 'n': case 'N':
            params[P_FIRSTDEEP] = readInt(optName, optArg);
            if (params[P_FIRSTDEEP] <= 0)
               printError("first depth must be positive.");
            break;
         case 'c': case 'C':
            params[P_CACHEMEM] = readInt(optName, optArg);
            break;
         case 'i': case 'I':
            params[P_MINDEEP] = readInt(optName, optArg);
            break;
         case 'q': case 'Q':
            params[P_QBITS] = readInt(optName, optArg);
            break;
         case 'h': case 'H':
            params[P_HASHBITS] = readInt(optName, optArg);
            break;
         case 'b': case 'B':
            params[P_BASEBITS] = readInt(optName, optArg);
            break;
#ifdef _OPENMP
         case 't': case 'T':
            params[P_NUMTHREADS] = readInt(optName, optArg);
            break;
#endif
         case 'f': case 'F':
            params[P_NUMSHIPS] = readInt(optName, optArg);
            break;
         case 'g': case 'G':
            params[P_MINEXTENSION] = readInt(optName, optArg);
            break;
         case 'p': case 'P':
            previewFlag = 1;
            break;
         case 'd': case 'D':
            dumpRoot = optArg;
            if (strlen(dumpRoot) > MAXDUMPROOT)
               printError("dump root exceeds maximum allowed length (" XSTR(MAXDUMPROOT) ")");
            break;
         case 'j': case 'J':
            splitNum = readInt(optName, optArg);
            if (splitNum < 0) splitNum = 0;
            break;
         case 'e': case 'E':
            initRows = optArg;
            initRowsFlag = 1;
            break;
         case 'l': case 'L':
            loadFile = optArg;
            loadDumpFlag = 1;
            loadParams();
            break;
         case 'a': case 'A':
            params[P_DUMPINTERVAL] = readInt(optName, optArg);
            if (params[P_DUMPINTERVAL] < 0)
               printError("dump interval must be nonnegative");
            break;
         case 259:   /* --enable-deep-printing */
            params[P_PRINTDEEP] = 1;
            break;
         case 260:   /* --disable-deep-printing */
            params[P_PRINTDEEP] = 0;
            break;
         case 261:   /* --enable-longest-partial */
            params[P_LONGEST] = 1;
            break;
         case 262:   /* --disable-longest-partial */
            params[P_LONGEST] = 0;
            break;
         case 263:   /* --dump-mode */
            switch(optArg[0]) {
               case 'o': case 'O':
                  params[P_DUMPMODE] = D_OVERWRITE; break;
               case 's': case 'S':
                  params[P_DUMPMODE] = D_SEQUENTIAL; break;
               case 'd': case 'D':
                  params[P_DUMPMODE] = D_DISABLED; break;
               default:
                  optError("unrecognized dump mode ", optArg);
                  break;
            }
            break;
         case 264:   /* --fixed-depth */
            params[P_EVERYDEPTH] = 1;
            params[P_MINDEEP] = 1;
            params[P_FIRSTDEEP] = readInt(optName, optArg);
            if (params[P_FIRSTDEEP] <= 0)
               printError("fixed depth must be positive.");
            break;
         case 265:   /* --enable-early-exit */
            params[P_EARLYEXIT] = 1;
            break;
         case 268:   /* --tables auto|pair|factored */
            if (!optArg) printError("missing argument for --tables");
            else if (optArg[0] == 'a') tableMode = 0;
            else if (optArg[0] == 'p') tableMode = 1;
            else if (optArg[0] == 'f') tableMode = 2;
            else optError("unrecognized table mode ", optArg);
            break;
         case 270:   /* --lookahead-budget N */
            lookaheadBudget = readInt(optName, optArg);
            if (lookaheadBudget < 0) printError("lookahead budget must be nonnegative.");
            break;
         case 271:   /* --split-cache-bits N */
            splitCacheBits = readInt(optName, optArg);
            if (splitCacheBits < 8 || splitCacheBits > 30) printError("split cache bits must be between 8 and 30.");
            break;
         case 269:   /* --lookahead-method auto|cone|split|hybrid */
            if (!optArg) printError("missing argument for --lookahead-method");
            else if (optArg[0] == 'a') lookaheadMethod = 0;
            else if (optArg[0] == 'c') lookaheadMethod = 1;
            else if (optArg[0] == 's') lookaheadMethod = 2;
            else if (optArg[0] == 'h') lookaheadMethod = 3;
            else optError("unrecognized lookahead method ", optArg);
            break;
         case 267:   /* --lookahead-depth */
            lookaheadDepth = readInt(optName, optArg);
            if (lookaheadDepth < 1 || lookaheadDepth > 3)
               printError("lookahead depth must be 1, 2 or 3.");
            break;
         case 266:   /* --disable-early-exit */
            params[P_EARLYEXIT] = 0;
            break;
         case 256:   /* --help */
            printHelp();
            break;
         case ':':
            optError("missing argument for option ", optName);
            break;
         case '?':
            optError("unrecognized option ", optName);
            break;
         default:
            printError("option parser failed");
            break;
      }
   }
}

/* ========================================= */
/*  Set up search with the given parameters  */
/* ========================================= */

void searchSetup(void) {
   if (params[P_CACHEMEM] < 0){
      if (5 * params[P_OFFSET] > params[P_PERIOD]) params[P_CACHEMEM] *= -1;
      else params[P_CACHEMEM] = 0;
   }
   
   checkParams();  /* Exit if parameters are invalid */
   
   if (aborting){
      fprintf(stderr, "\nUse --help for a list of available options.\n");
      exit(1);
   }
   
   if (loadDumpFlag) loadState();
   else {
      width = params[P_WIDTH];
      period = params[P_PERIOD];
      offset = params[P_OFFSET];
      
      nRowsInState = period+period;
      
      params[P_DEPTHLIMIT] = DEFAULT_DEPTHLIMIT;
      
      queueBase = (node*)malloc((QSIZE>>BASEBITS)*sizeof(node));
      queueRows = (row*)malloc(QSIZE*sizeof(row));
      if (queueBase == 0 || queueRows == 0) {
         printf("Unable to allocate BFS queue!\n");
         exit(1);
      }
      
      if (hashBits == 0) seenHash = 0;
      else {
         seenHash = (node*)malloc(HASHSIZE*sizeof(node));
         if (seenHash == 0) printf("Unable to allocate hash table, duplicate elimination disabled\n");
      }
      
      extensionRows = (row**)calloc(1LLU << (params[P_DEPTHLIMIT] + 1), sizeof(*extensionRows));
      extensionIndex = (uint32_t*)calloc(EXTINDEX_SIZE, sizeof(*extensionIndex));
      
      resetQ();
      resetHash();
      
      enqueue(0,0);
      
      if (initRowsFlag) loadInitRows(initRows);
   }
   
#ifndef QSIMPLE
   makePhases();
   makeSubperiodTables();
#endif
   
   /* Generate proper rule string for printing patterns */
   uint32_t i;
   uint32_t j = 0;
   int k = 1;
   for (i = 0; i < 151 && rule[i] != '\0'; i++){
      if (rule[i] == '~') k = 0;
      else if (rule[i] == '/') k = 1;
      if (k) baseRule[j++] = rule[i];
   }
   baseRule[j] = '\0';
   
   if (previewFlag){
      params[P_NUMSHIPS] = 0;
      preview();
      exit(0);
   }
   
   if (params[P_MINEXTENSION] < 0) params[P_MINEXTENSION] = 0;
   if (params[P_FIRSTDEEP] < 0) params[P_FIRSTDEEP] = 0;
   
   dumpMode = params[P_DUMPMODE];   /* these may need to be different when splitting */
   
   /* split queue across multiple files */
   if (splitNum > 0){
      dumpMode = D_SEQUENTIAL; 
      
      echoParams();
      printf("\n");
      
      if (!loadDumpFlag || qHead == 0 || splitNum == 1){
         dumpFlag = DUMPPENDING;
         if (qHead == 0){      /* can't use doCompact() here, because it tries to access rows[-1] */
            qStart = qHead;
            qEnd = qTail;
            dumpState();
         }
         else doCompact();
         if (dumpFlag == DUMPSUCCESS){
            printf("State dumped to %s\n",dumpFile);
            exit(0);
         }
         else{
            fprintf(stderr, "Error: dump failed.\n");
            exit(1);
         }
      }
      
      if (strlen(dumpRoot) > 244){
         fprintf(stderr, "Error: dump root too long.\n");
         exit(1);
      }
      
      if (splitNum >= 100000){
         fprintf(stderr, "Warning: queue cannot be split into more than 99999 files.\n");
         splitNum = 99999;
      }
      
      /* count nodes in queue */
      int totalNodes = 0;
      node x;
      for (x = qHead; x < qTail; x++){
         if (!EMPTY(x)) totalNodes++;
      }
      
      /* nodes per file is rounded up */
      unsigned long nodesPerFile = (totalNodes - 1) / splitNum + 1;
      
      printf("Splitting search state with %lu queue nodes per file\n",nodesPerFile);
      
      /* save qHead and qTail, as creating the pieces will change their values */
      node fixedQHead = qHead;
      node fixedQTail = qTail;
      
      /* delete the queue; we will reload it as needed */
      free(queueBase);
      free(queueRows);
      free(seenHash);
      
      uint32_t deepIndex;
      for (deepIndex = 0; deepIndex < 1LLU << (params[P_DEPTHLIMIT] + 1); ++deepIndex){
         if (extensionRows[deepIndex]) free(extensionRows[deepIndex]);
         extensionRows[deepIndex] = 0;
      }
      free(extensionRows);
      free(extensionIndex);
      
      node currNode = fixedQHead;
      
      dumpNum = 1;   /* reset dumpNum prior to splitting */
      int firstDumpNum = 0;
      while (currNode < fixedQTail){
         
         /* load the queue */
         loadState();
         
         /* empty everything before currNode */
         j = extQHead;
         for (x = fixedQHead; x < currNode; ++x){
            MAKEEMPTY(x);
            extensionIndex[j] = 0;
            ++j;
         }
         
         /* retain the specified number of nonempty nodes */
         i = 0;
         while (i < nodesPerFile && x < fixedQTail){
            if (!EMPTY(x))
               ++i;
            ++x;
            ++j;
         }
         
         /* update currNode */
         currNode = x;
         
         /* empty everything after specified nonempty nodes */
         while (x < fixedQTail){
            MAKEEMPTY(x);
            extensionIndex[j] = 0;
            ++x;
            ++j;
         }

         /* save the piece */
         dumpFlag = DUMPPENDING;
         doCompact();
         
         if (!firstDumpNum) firstDumpNum = dumpNum - 1;
         
         if (dumpFlag != DUMPSUCCESS){
            printf("Failed to save %s\n",dumpFile);
            exit(1);
         }
         
         if (dumpNum >= DUMPLIMIT){
            fprintf(stderr, "Error: dump file number limit (" XSTR(DUMPLIMIT) ") reached.\n");
            fprintf(stderr, "       Try splitting the queue in a new directory.\n");
            exit(1);
         }
         
         /* free memory allocated in loadState() */
         for (deepIndex = 0; deepIndex < 1LLU << (params[P_DEPTHLIMIT] + 1); ++deepIndex){
            if (extensionRows[deepIndex]) free(extensionRows[deepIndex]);
            extensionRows[deepIndex] = 0;
         }
         
         free(queueBase);
         free(queueRows);
         free(seenHash);
         free(extensionRows);
         free(extensionIndex);
      }
      
      printf("Saved pieces in files %s%05d to %s\n",dumpRoot,firstDumpNum,dumpFile);
      exit(0);
   }
   
   omp_set_num_threads(params[P_NUMTHREADS]);
   
   memlimit = params[P_MEMLIMIT] >= 0 ? ((long long)params[P_MEMLIMIT]) << 20 : 0;
   
   /* Allocate lookahead cache */
#ifndef NOCACHE
   cachesize = 32768;
   while (cachesize * ((long long) sizeof(cacheentry)) < 550000 * (long long) params[P_CACHEMEM])
      cachesize <<= 1;
   memusage += sizeof(cacheentry) * (cachesize + 5) * params[P_NUMTHREADS];
   if (params[P_MEMLIMIT] >= 0 && memusage > memlimit){
      printf("Not enough memory to allocate lookahead cache\n");
      exit(1);
   }
   totalCache = (cacheentry *)calloc((size_t) ((cachesize + 5) * params[P_NUMTHREADS]), 
                                     sizeof(cacheentry));
   cache = (cacheentry **)calloc(params[P_NUMTHREADS], sizeof(*cache));
   
   for (k = 0; k < params[P_NUMTHREADS]; k++)
      cache[k] = totalCache + (cachesize + 5) * k;
#endif
   
   chooseTableMode();
   echoParams();
   
   makeRuleTable3x3();
   makeTables();
   
   rephase();
   
   parseDumpRoot();
   time(&lastDumpTime);
   
   timeStamp();
}

void finalReport(void) {
   timeStamp();
   printf("Search complete.\n\n");
   printf("Total time %.1f s, peak memory %.0f MB, tables+index %lld MB\n", elapsedSeconds(), peakRssMB(), memusage >> 20);
   
   printf("%d %s%s found.\n",numFound,(params[P_BOUNDARYSYM] == SYM_UNDEF) ? "spaceship" : "wave",(numFound == 1) ? "" : "s");
   printf("Maximum depth reached: %d\n",longest);
   if (params[P_LONGEST] && aborting != 3){ /* aborting == 3 means we reached ship limit */
      if (patternBuf) printf("Longest partial result:\n\n%s",patternBuf);
      else printf("No partial results found.\n");
   }
}
