/* qfind --- A spaceship search program by Matthias Merzenich.
** 
** This is an attempt at combining the functionality of gfind and zfind.
**
** Based on code by David Eppstein, zdr, Paul Tooke, and Tomas Rokicki.
** Thanks also to Frank Everdij, Adam P. Goucher, Alex Greason, and
** praosylen for additional code and suggestions.
**
** How the search works (overview added with the optimized version):
**  - A partial spaceship is a sequence of rows.  Because the ship moves `offset`
**    cells every `period` generations, consecutive rows of the sequence are rows of
**    different generations, and each new row is constrained by the rule applied to
**    rows 1 and 2 periods earlier.  Lookup tables (pairTable, built per row pair on
**    demand by buildTable() in common.h) list, for rows (row1,row2) and a required
**    result row, every row3 that produces it.  getSuccessors() reads such a list.
**  - Breadth-first: process() takes a node from the queue and enqueues each candidate
**    next row that passes lookAhead() and was not seen before (seenHash).
**  - When the queue fills up, deepen() runs depthFirst() on every queue node in
**    parallel: a node survives only if it can be extended by the current deepening
**    amount; dead nodes are removed and the queue is compacted (doCompact()).
**  - lookAhead() is a cheap necessary test that rejects rows that cannot be extended
**    (see the block above lookAheadReference()).
**  - A ship is reported when a branch ends in enough empty rows (terminal()).
*/

/*  qfind allows the spaceship velocity to be set at run time.  However,     */
/*  the option is available to compile qfind with a predetermined velocity.  */
/*  This may result in the search running slightly faster.  To compile with  */
/*  a predetermined velocity, uncomment the two lines below and change the   */
/*  values as desired.  For searches with a predetermined velocity you must  */
/*  have gcd(PERIOD,OFFSET) = 1.                                             */

//#define PERIOD 5
//#define OFFSET 2

/*  Lookahead caching seems to speed up searches for spaceships with         */
/*  velocities exceeding c/5 and slow down searches for spaceships with      */
/*  velocities at or below c/5.  If you compile qfind with a predetermined   */
/*  velocity, lookahead caching will automatically be enabled or disabled    */
/*  depending on your choice of PERIOD and OFFSET.  To override the          */
/*  automatic choice, uncomment one of the following two lines.              */

//#define NOCACHE
//#define FORCECACHE

#if !defined(PERIOD) && !defined(OFFSET)
   #define PERIOD period
   #define OFFSET offset
   #define FWDOFF(x) fwdOff[x]
   #define BACKOFF(x) backOff[x]
   #define DOUBLEOFF(x) doubleOff[x]
   #define TRIPLEOFF(x) tripleOff[x]
#elif defined(PERIOD) && !defined(OFFSET)
   #error "OFFSET must be defined."
   #define OFFSET period   /* prevents irrelevant compiler errors */
#elif !defined(PERIOD) && defined(OFFSET)
   #error "PERIOD must be defined."
   #define PERIOD period   /* prevents irrelevant compiler errors */
#else
   #define QSIMPLE

   #define FWDOFF(x) (OFFSET)
   #define BACKOFF(x) (OFFSET)
   #define DOUBLEOFF(x) (2*OFFSET)
   #define TRIPLEOFF(x) (3*OFFSET)

   /* If QSIMPLE is defined, the last argument of lookAhead() is unused.  */
   /* The following line prevents an "unused parameter" compiler warning. */
   #define lookAhead(x,y,z) lookAhead(x,y)
   
   #if PERIOD < 1 || OFFSET < 1
      #error "Invalid value for PERIOD or OFFSET."
   #endif
   
   #if OFFSET > PERIOD && PERIOD > 0
      #error "OFFSET cannot exceed PERIOD."
   #endif
   
   #if OFFSET == PERIOD && PERIOD > 0
      #error "Photons are not supported."
   #endif
   
   #if 2 * OFFSET > PERIOD && PERIOD > 0
      #warning "Searches for speeds exceeding c/2 may not work correctly."
   #endif
   
   /* The uncached fast lookahead measured faster than the lookahead cache at every */
   /* speed tested, so the cache is compiled out unless FORCECACHE is defined.     */
   #ifndef FORCECACHE
      #define NOCACHE
   #endif
#endif

/* FWDOFF is undefined only if there is already an error in compilation. */
/* The following macro definitions suppress unwanted compiler warnings.  */
#ifndef FWDOFF
   #define FWDOFF(x) rowPhase   /* Prevents unused parameter warning */
   #define BACKOFF(x) 0
   #define DOUBLEOFF(x) 0
   #define TRIPLEOFF(x) 0
#endif

#include "common.h"

#ifdef QSIMPLE
   #define ROWPHASE_ARG 0
#else
   #define ROWPHASE_ARG rowPhase
#endif

#ifdef QSIMPLE
   #if PERIOD > MAXPERIOD
      #error "maximum allowed PERIOD exceeded."
   #endif
#endif

#ifndef QSIMPLE
/* Phase-dependent index offsets in the row sequence (computed by makePhases()):
**   backOff[ph]:  a new row at index a (phase ph) is chosen so that the rows at
**                 a-2P, a-P and a (stacked) evolve the middle row a-P into row
**                 a-P+backOff[ph]; this is the lookup getSuccessors() performs.
**   fwdOff[ph], doubleOff[ph], tripleOff[ph]: distances to the rows 1, 2 and 3
**                 generation steps away, used by the lookahead (see the notation
**                 block above lookAheadReference()).
** With a fixed speed (QSIMPLE) these are the constants OFFSET, 2*OFFSET, 3*OFFSET. */
int fwdOff[MAXPERIOD], backOff[MAXPERIOD], doubleOff[MAXPERIOD], tripleOff[MAXPERIOD];

void makePhases(void){
   int i;
   for (i = 0; i < period; i++) backOff[i] = -1;
   i = 0;
   for (;;) {
      int j = offset;
      while (backOff[(i+j)%period] >= 0 && j < period) j++;
      if (j == period) {
         backOff[i] = period-i;
         break;
      }
      backOff[i] = j;
      i = (i+j)%period;
   }
   for (i = 0; i < period; i++)
      fwdOff[(i+backOff[i])%period] = backOff[i];
   for (i = 0; i < period; i++) {
      int j = (i - fwdOff[i]);
      if (j < 0) j += period;
      doubleOff[i] = fwdOff[i] + fwdOff[j];
   }
   for (i = 0; i <  period; i++){
      int j = (i - fwdOff[i]);
      if (j < 0) j += period;
      tripleOff[i] = fwdOff[i] + doubleOff[j];
   }
}
#endif

/* ========================================================================== *
** Lookahead
**
** lookAhead(pathRows, a, rowPhase) decides whether the candidate row pathRows[a]
** can possibly be extended.  Notation: "rowLG" is the row L periods further
** along the ship than the candidate's period block, G generations ahead of it;
** listLG is the list of values rowLG may take, read from the lookup tables.
** The candidate is kept iff there exist
**
**      row11 in list11, row12 in list12, row13 in list13,
**      row22 in list22(row12,row11), row23 in list23(row13,row12)
**
** such that the table for (row13,row23) has a non-empty bucket for row22, i.e.
** some row33 exists.  This 3-generation "cone" is a necessary condition for
** any extension, so dropping candidates that fail it never loses a ship.
**
** The rows the tables are looked up by (all at indices below a, i.e. already
** fixed) are:
**      rowA1 = pathRows[a-P-f]  rowA2 = pathRows[a-f]     (generation 1)
**      rowB1 = pathRows[a-P-d]  rowB2 = pathRows[a-d]     (generation 2)
**      rowC1 = pathRows[a-P-t]  rowC2 = pathRows[a-t]     (generation 3)
** with P = PERIOD and f, d, t = FWDOFF, DOUBLEOFF, TRIPLEOFF of rowPhase.
** If t >= P, row13 is itself an already fixed row (pathRows[a+P-t]); then
** list13 has exactly that one entry ("row13Known").
** ========================================================================== */

/* Decode buffers (one list each) used by lookAheadReference(). */
static _Thread_local uint16_t *refDecodeBuf[5];
static inline uint16_t *decodeBuffer(uint16_t **slot) {
   if (__builtin_expect(!*slot, 0)) *slot = (uint16_t*)malloc(DECODE_CAP * sizeof(uint16_t));
   return *slot;
}

/* The original lookahead, unchanged in logic.  Used when the lookahead cache is
** enabled (speeds above c/5 by default), when TRIPLEOFF == PERIOD, and as the
** reference that -DLA_VERIFY compares lookAheadFast() against.               */
__attribute__((unused))
static int lookAheadReference(row *pathRows, int candRow, int rowPhase){
   (void)rowPhase;
   int i11, i12, i13, i22, i23;
   const uint16_t *list11, *list12, *list13, *list22, *list23;
   uint64_t id11, id12, id13 = 0;     /* list identities, used as cache keys */
   (void)id11; (void)id12; (void)id13;
   int len11, len12, len13, len22, len23;
   int row11, row12, row13, row22, row23;
#ifndef NOCACHE
   int cacheSlot;
#endif
   
   getSuccessorsWithId(pathRows[candRow - PERIOD - FWDOFF(rowPhase)],
                       pathRows[candRow - FWDOFF(rowPhase)],
                       pathRows[candRow], &list11, &len11, decodeBuffer(&refDecodeBuf[0]), &id11);
   if (!len11)
      return 0;
   
   getSuccessorsWithId(pathRows[candRow - PERIOD - DOUBLEOFF(rowPhase)],
                       pathRows[candRow - DOUBLEOFF(rowPhase)],
                       pathRows[candRow - FWDOFF(rowPhase)], &list12, &len12, decodeBuffer(&refDecodeBuf[1]), &id12);
   /* An empty list12 or list13 means "no".  Return before using the lookahead cache: */
   /* an empty list has no unique identity (in the original code its pointer equals   */
   /* the start of the next non-empty list), so it must never become part of a key.  */
   if (!len12)
      return 0;
   
#if !defined(QSIMPLE)   /* No need for conditional when using QSIMPLE */
   if (TRIPLEOFF(rowPhase) >= PERIOD)
#endif
   
#if !defined(QSIMPLE) || 3*OFFSET >= PERIOD
   {
      list13 = pathRows + (candRow + PERIOD - TRIPLEOFF(rowPhase));
      len13 = 1;
   #ifndef NOCACHE
      id13 = (uint64_t)(uintptr_t)(rowOrderWeight + list13[0]);   /* any unique non-zero value per row */
      cacheSlot = cacheLookup(id11, id12, id13,
         (pathRows[candRow-DOUBLEOFF(rowPhase)] << width) + pathRows[candRow-TRIPLEOFF(rowPhase)]);
   #endif
   }
#endif
   
#if !defined(QSIMPLE)   /* No need for conditional when using QSIMPLE */
   else
#endif
   
#if !defined(QSIMPLE) || 3*OFFSET < PERIOD
   {
      getSuccessorsWithId(pathRows[candRow - PERIOD - TRIPLEOFF(rowPhase)],
                          pathRows[candRow - TRIPLEOFF(rowPhase)],
                          pathRows[candRow - DOUBLEOFF(rowPhase)], &list13, &len13, decodeBuffer(&refDecodeBuf[2]), &id13);
      if (!len13)
         return 0;
   #ifndef NOCACHE
      cacheSlot = cacheLookup(id11, id12, id13,
         (pathRows[candRow-DOUBLEOFF(rowPhase)] << width) + pathRows[candRow-TRIPLEOFF(rowPhase)]);
   #endif
   }
#endif
   
#ifndef NOCACHE
   if (cacheSlot < 0)          /* cache hit: -2 = "fails", -1 = "passes" */
      return cacheSlot + 2;
#endif
   
   for (i11 = 0; i11 < len11; ++i11){
      row11 = list11[i11];
      for (i12 = 0; i12 < len12; ++i12){
         row12 = list12[i12];
         getSuccessors(pathRows[candRow - DOUBLEOFF(rowPhase)],
                       row12, row11, &list22, &len22, decodeBuffer(&refDecodeBuf[3]));
         if (!len22) continue;
         
         for (i13 = 0; i13 < len13; ++i13){
            row13 = list13[i13];
            getSuccessors(pathRows[candRow - TRIPLEOFF(rowPhase)],
                          row13, row12, &list23, &len23, decodeBuffer(&refDecodeBuf[4]));
            if (!len23) continue;
            
            for (i23 = 0; i23 < len23; ++i23){
               row23 = list23[i23];
               TableRef table33 = TABLE_REF(useFactored, row13, row23);
               for (i22 = 0; i22 < len22; ++i22){
                  row22 = list22[i22];
                  if (TABLE_HAS_ROW(useFactored, table33, row22)) {
#ifndef NOCACHE
                     cacheStore(cacheSlot, 1);
#endif
                     return 1;
                  }
               }
            }
         }
      }
   }
#ifndef NOCACHE
   cacheStore(cacheSlot, 0);
#endif
   return 0;
}


/* -------------------------------------------------------------------------- *
** Fast lookahead (used whenever the lookahead cache is disabled).
**
** Only list11 depends on the candidate row; list12, list13 and the tables behind
** list22/list23 depend only on rows above it.  The depth-first search tries many
** sibling candidates in a row, so those parts are kept per search-tree level in
** a LookaheadLevel and reused by every sibling:
**   stage 0 (laBuildStage0): the looked-up rows and the table giving list11.
**        About 70% of candidates are rejected here by one bit test (no row11).
**   stage 1 (laBuildStage1): list12, list13, and caches for the per-row12 and
**        per-row13 tables (filled lazily).
** The existence test itself is the same as lookAheadReference(), only the loop
** order differs (row12 outermost), so the result is identical for every input.
** A level's data is discarded (laForgetLevel) whenever the search enters that
** level with a new parent, and for all levels (laForgetAll) at the start of each
** depth-first search.
** -------------------------------------------------------------------------- */
typedef struct {
   int built0, built1;           /* stage 0 / stage 1 data valid                    */
   int rowA1, rowA2, rowB1, rowB2, rowC1, rowC2;   /* see notation above            */
   int row13Known;               /* TRIPLEOFF >= PERIOD: row13 is the fixed row rowC1 */
   TableRef table11;             /* table whose bucket [candidate] is list11        */
   const uint16_t *list12; int len12;
   const uint16_t *list13; int len13;
   uint16_t only13;              /* list13 storage when row13Known                  */
   uint16_t *decode12, *decode13;   /* decode buffers for list12/list13 (-DPACKLISTS) */
   TableRef *tables22;           /* table(rowB2,row12) per entry of list12, 0 = not yet fetched */
   TableRef *tables23;           /* table(rowC2,row13) per entry of list13, 0 = not yet fetched */
   int cap12, cap13;             /* allocated sizes of tables22 / tables23         */
} LookaheadLevel;

static _Thread_local LookaheadLevel **laLevels;
static _Thread_local int laLevelCap;
static _Thread_local const uint16_t **laLists22;   /* list22(row12,row11) for each row11 */
static _Thread_local int *laLens22;
static _Thread_local int laLists22Cap;
static _Thread_local uint16_t *laDecode11, *laDecode23, *laDecode22;   /* decode buffers */

static void laGrowLevels(int level) {
   int ncap = laLevelCap ? laLevelCap : 64;
   while (ncap <= level) ncap *= 2;
   laLevels = (LookaheadLevel**)realloc(laLevels, ncap * sizeof(*laLevels));
   for (int i = laLevelCap; i < ncap; ++i) laLevels[i] = (LookaheadLevel*)calloc(1, sizeof(LookaheadLevel));
   laLevelCap = ncap;
}

static inline LookaheadLevel *laGetLevel(int level) {
   if (__builtin_expect(level >= laLevelCap, 0)) laGrowLevels(level);
   return laLevels[level];
}

/* forget the data of one level / of all levels */
static void splitForgetLevel(int level);
static void splitForgetAll(void);
static inline void laForgetLevel(int level) {
   LookaheadLevel *lvl = laGetLevel(level);
   lvl->built0 = lvl->built1 = 0;
   splitForgetLevel(level);
}
static void laForgetAll(void) {
   for (int i = 0; i < laLevelCap; ++i) laLevels[i]->built0 = laLevels[i]->built1 = 0;
   splitForgetAll();
}

static void laBuildStage0(LookaheadLevel *lvl, const row *pathRows, int candRow, int rowPhase) {
   (void)rowPhase;
   const int f = FWDOFF(rowPhase), d = DOUBLEOFF(rowPhase), t = TRIPLEOFF(rowPhase);
   lvl->row13Known = (t >= PERIOD);
   lvl->rowA1 = pathRows[candRow - PERIOD - f]; lvl->rowA2 = pathRows[candRow - f];
   lvl->rowB1 = pathRows[candRow - PERIOD - d]; lvl->rowB2 = pathRows[candRow - d];
   lvl->rowC2 = pathRows[candRow - t];
   lvl->rowC1 = lvl->row13Known ? pathRows[candRow + PERIOD - t] : pathRows[candRow - PERIOD - t];
   lvl->table11 = TABLE_REF(useFactored, lvl->rowA1, lvl->rowA2);
   lvl->built0 = 1;
}

static void laBuildStage1(LookaheadLevel *lvl) {
   if (!lvl->decode12) {
      lvl->decode12 = (uint16_t*)malloc(DECODE_CAP * sizeof(uint16_t));
      lvl->decode13 = (uint16_t*)malloc(DECODE_CAP * sizeof(uint16_t));
   }
   lvl->list12 = TABLE_ROWS(useFactored, TABLE_REF(useFactored, lvl->rowB1, lvl->rowB2), lvl->rowA2, &lvl->len12, lvl->decode12);
   if (lvl->row13Known) { lvl->only13 = (uint16_t)lvl->rowC1; lvl->list13 = &lvl->only13; lvl->len13 = 1; }
   else lvl->list13 = TABLE_ROWS(useFactored, TABLE_REF(useFactored, lvl->rowC1, lvl->rowC2), lvl->rowB2, &lvl->len13, lvl->decode13);
   if (lvl->len12 > lvl->cap12) {
      lvl->cap12 = 2 * lvl->len12;
      lvl->tables22 = (TableRef*)realloc((void*)lvl->tables22, lvl->cap12 * sizeof(*lvl->tables22));
   }
   if (lvl->len13 > lvl->cap13) {
      lvl->cap13 = 2 * lvl->len13;
      lvl->tables23 = (TableRef*)realloc((void*)lvl->tables23, lvl->cap13 * sizeof(*lvl->tables23));
   }
   if (lvl->len12) memset(lvl->tables22, 0, lvl->len12 * sizeof(*lvl->tables22));
   if (lvl->len13) memset(lvl->tables23, 0, lvl->len13 * sizeof(*lvl->tables23));
   lvl->built1 = 1;
}

/* Compiled twice (factored = 0: per-pair tables, 1: factored tables), see lookAheadFast() */
/* Work budget for the hybrid lookahead: lookAheadFastT() gives up (returns -1)
** once it has done more than laBudget units of work (list entries scanned and row33
** tests); lookAhead() then answers the call with lookAheadSplit(), whose cost does
** not grow with the width.  Both are exact, so the answer is the same either way.
** laBudget < 0 means unlimited.                                                    */
static _Thread_local long laBudget = -1;


#ifdef INNER_SPLIT2
/* Like factoredInner(), but row11 is also split: its low half only affects the low
** side and its high half only the high side, so within each row11 overlap value w the
** low masks are ORed over all low halves of row11 and the high masks over all high
** halves (cost |L11|+|H11| instead of |L11|*|H11|).  table11 = factored ref of
** (rowA1,rowA2), cand = the candidate row (row11's result).                       */
static void innerAcc(int isHigh, const uint8_t *e11, int s11, int t11, uint32_t p22, uint32_t p23,
                     int r12, int x13, int m23, uint8_t acc[16]) {
   const int lb = lowHalfBits, hb = highHalfBits, lr = lowResultBits, hr = highResultBits, rm = (1 << lr) - 1;
   const int bits = isHigh ? hb : lb, res = isHigh ? hr : lr;
   const uint8_t *masks = isHigh ? highHalfMasks : lowHalfMasks, *rows = isHigh ? highHalfRows : lowHalfRows;
   const uint16_t *offs = isHigh ? highHalfOffsets : lowHalfOffsets;
   const size_t opp = isHigh ? HIGH_OFFSETS_PER_PAIR : LOW_OFFSETS_PER_PAIR;
   const uint16_t *o23 = offs + (size_t)p23 * opp + (r12 << 2);
   const uint8_t *e23 = rows + ((size_t)p23 << bits), *e22 = rows + ((size_t)p22 << bits);
   for (int i11 = s11; i11 < t11; i11++) {
      const int r11 = isHigh ? (e11[i11] >> 1) : (e11[i11] & rm);
      const int m22 = masks[((size_t)p22 << res) | r11];
      if (!m22) continue;
      const uint16_t *o22 = offs + (size_t)p22 * opp + (r11 << 2);
      for (int u = 0; u < 4; u++) {
         if (!((m23 >> u) & 1)) continue;
         for (int i = o23[u]; i < o23[u + 1]; i++) {
            const uint8_t *mk = masks + ((size_t)((x13 << bits) | e23[i]) << res);
            for (int v = 0; v < 4; v++) {
               if (!((m22 >> v) & 1)) continue;
               int a = 0;
               for (int j = o22[v]; j < o22[v + 1]; j++) a |= mk[isHigh ? (e22[j] >> 1) : (e22[j] & rm)];
               acc[(u << 2) | v] |= (uint8_t)a;
            }
         }
      }
   }
}
static int factoredInner2(uint32_t ref11, int cand, uint32_t ref22, uint32_t ref23, int row12, int row13) {
   const int lb = lowHalfBits, hb = highHalfBits, hs = highHalfShift, lr = lowResultBits, hr = highResultBits, rm = (1 << lr) - 1;
   const uint32_t p11L = ref11 & 0xFFFF, p11H = ref11 >> 16, p22L = ref22 & 0xFFFF, p22H = ref22 >> 16, p23L = ref23 & 0xFFFF, p23H = ref23 >> 16;
   const int cL = cand & rm, cH = cand >> lr, r12L = row12 & rm, r12H = row12 >> lr;
   const int m11 = lowHalfMasks[((size_t)p11L << lr) | cL] & highHalfMasks[((size_t)p11H << hr) | cH];
   const int m23 = lowHalfMasks[((size_t)p23L << lr) | r12L] & highHalfMasks[((size_t)p23H << hr) | r12H];
   if (!m11 || !m23) return 0;
   const uint16_t *o11L = lowHalfOffsets + (size_t)p11L * LOW_OFFSETS_PER_PAIR + (cL << 2);
   const uint16_t *o11H = highHalfOffsets + (size_t)p11H * HIGH_OFFSETS_PER_PAIR + (cH << 2);
   const uint8_t *e11L = lowHalfRows + ((size_t)p11L << lb), *e11H = highHalfRows + ((size_t)p11H << hb);
   const int l13 = row13 & ((1 << lb) - 1), h13 = row13 >> hs;
   for (int w = 0; w < 4; w++) {
      if (!((m11 >> w) & 1)) continue;
      uint8_t aL[16] = {0}, aH[16] = {0};
      innerAcc(0, e11L, o11L[w], o11L[w + 1], p22L, p23L, r12L, l13, m23, aL);
      int any = 0; for (int k = 0; k < 16; k++) any |= aL[k];
      if (!any) continue;
      innerAcc(1, e11H, o11H[w], o11H[w + 1], p22H, p23H, r12H, h13, m23, aH);
      for (int k = 0; k < 16; k++) if (aL[k] & aH[k]) return 1;
   }
   return 0;
}
#endif
#ifdef INNER_SPLIT
#define DECODE_CAP_MAX ((1 << 14) + 8)
/* Innermost cone test in factored mode, on half-tables: does some row23 in
** list23(row13,row12) and row22 in list22(row12,row11) give a row33?  Both lists are
** joins (low x high, matched on overlap), so with 4-bit row33-overlap masks the test is
** "exists u,v: OR(low masks) & OR(high masks) != 0": cost |L23||L22|+|H23||H22|
** instead of |list23||list22|.                                                     */
static int factoredInner(uint32_t ref22, int row11, uint32_t ref23, int row12, int row13) {
   const int lb = lowHalfBits, hb = highHalfBits, hs = highHalfShift, lr = lowResultBits, hr = highResultBits, rm = (1 << lr) - 1;
   const uint32_t p22L = ref22 & 0xFFFF, p22H = ref22 >> 16, p23L = ref23 & 0xFFFF, p23H = ref23 >> 16;
   const int r11L = row11 & rm, r11H = row11 >> lr, r12L = row12 & rm, r12H = row12 >> lr;
   const int m22 = lowHalfMasks[((size_t)p22L << lr) | r11L] & highHalfMasks[((size_t)p22H << hr) | r11H];
   const int m23 = lowHalfMasks[((size_t)p23L << lr) | r12L] & highHalfMasks[((size_t)p23H << hr) | r12H];
   if (!m22 || !m23) return 0;
   const uint16_t *o22L = lowHalfOffsets + (size_t)p22L * LOW_OFFSETS_PER_PAIR + (r11L << 2);
   const uint16_t *o22H = highHalfOffsets + (size_t)p22H * HIGH_OFFSETS_PER_PAIR + (r11H << 2);
   const uint16_t *o23L = lowHalfOffsets + (size_t)p23L * LOW_OFFSETS_PER_PAIR + (r12L << 2);
   const uint16_t *o23H = highHalfOffsets + (size_t)p23H * HIGH_OFFSETS_PER_PAIR + (r12H << 2);
   const uint8_t *e22L = lowHalfRows + ((size_t)p22L << lb), *e22H = highHalfRows + ((size_t)p22H << hb);
   const uint8_t *e23L = lowHalfRows + ((size_t)p23L << lb), *e23H = highHalfRows + ((size_t)p23H << hb);
   const int l13 = row13 & ((1 << lb) - 1), h13 = row13 >> hs;
   for (int u = 0; u < 4; u++) {
      if (!((m23 >> u) & 1)) continue;
      for (int v = 0; v < 4; v++) {
         if (!((m22 >> v) & 1)) continue;
         int aL = 0, aH = 0;
         for (int i = o23L[u]; i < o23L[u + 1]; i++) {
            const uint8_t *mk = lowHalfMasks + ((size_t)((l13 << lb) | e23L[i]) << lr);
            for (int j = o22L[v]; j < o22L[v + 1]; j++) aL |= mk[e22L[j] & rm];
         }
         if (!aL) continue;
         for (int i = o23H[u]; i < o23H[u + 1]; i++) {
            const uint8_t *mk = highHalfMasks + ((size_t)((h13 << hb) | e23H[i]) << hr);
            for (int j = o22H[v]; j < o22H[v + 1]; j++) aH |= mk[e22H[j] >> 1];
         }
         if (aL & aH) return 1;
      }
   }
   return 0;
}
#endif
static inline __attribute__((always_inline)) int lookAheadFastT(row *pathRows, int candRow, int rowPhase, const int factored){
#ifndef QSIMPLE
   /* list13 would depend on the candidate itself; use the reference code for that case. */
   if (TRIPLEOFF(rowPhase) == PERIOD) return lookAheadReference(pathRows, candRow, rowPhase);
#elif 3*OFFSET == PERIOD
   return lookAheadReference(pathRows, candRow, rowPhase);
#endif
   LookaheadLevel *lvl = laGetLevel(candRow);
   if (__builtin_expect(!lvl->built0, 0)) laBuildStage0(lvl, pathRows, candRow, rowPhase);
   
   /* stage 0: does the candidate have any row11 at all? */
   TableRef table11 = lvl->table11;
   const int cand = pathRows[candRow];
   if (!TABLE_HAS_ROW(factored, table11, cand)) return 0;
   if (__builtin_expect(!laDecode11, 0)) {
      laDecode11 = (uint16_t*)malloc(DECODE_CAP * sizeof(uint16_t));
      laDecode23 = (uint16_t*)malloc(DECODE_CAP * sizeof(uint16_t));
      laDecode22 = (uint16_t*)malloc(DECODE_CAP * sizeof(uint16_t));
   }
   int len11;
   const uint16_t *list11 = TABLE_ROWS(factored, table11, cand, &len11, laDecode11);
   
   if (!lvl->built1) laBuildStage1(lvl);
   
   if (len11 > laLists22Cap) {
      laLists22Cap = 2 * len11;
      laLists22 = (const uint16_t**)realloc((void*)laLists22, laLists22Cap * sizeof(*laLists22));
      laLens22 = (int*)realloc(laLens22, laLists22Cap * sizeof(*laLens22));
   }
   const uint16_t **lists22 = laLists22;
   int *lens22 = laLens22;
   long work = 0;
   const long budget = laBudget;
#ifdef INNER_SPLIT
   if (factored && lookaheadDepth == 3 && budget < 0) {
      for (int i12 = 0; i12 < lvl->len12; ++i12) {
         const int row12 = lvl->list12[i12];
         TableRef table22 = lvl->tables22[i12];
         if (!table22) table22 = lvl->tables22[i12] = TABLE_REF(factored, lvl->rowB2, row12);
         uint16_t keep11[DECODE_CAP_MAX]; int nk = 0;   /* row11s with a non-empty list22 */
         for (int k = 0; k < len11; ++k)
            if (FACTORED_HAS_ROW(table22, list11[k])) keep11[nk++] = list11[k];
         if (!nk) continue;
         for (int i13 = 0; i13 < lvl->len13; ++i13) {
            const int row13 = lvl->list13[i13];
            TableRef table23 = lvl->tables23[i13];
            if (!table23) table23 = lvl->tables23[i13] = TABLE_REF(factored, lvl->rowC2, row13);
            if (!FACTORED_HAS_ROW(table23, row12)) continue;
#ifdef INNER_SPLIT2
            (void)nk;
            if (factoredInner2((uint32_t)table11, cand, (uint32_t)table22, (uint32_t)table23, row12, row13)) return 1;
#else
            for (int k = 0; k < nk; ++k)
               if (factoredInner((uint32_t)table22, keep11[k], (uint32_t)table23, row12, row13)) return 1;
#endif
         }
      }
      return 0;
   }
#endif
   
   for (int i12 = 0; i12 < lvl->len12; ++i12) {
      const int row12 = lvl->list12[i12];
      TableRef table22 = lvl->tables22[i12];
      if (!table22) table22 = lvl->tables22[i12] = TABLE_REF(factored, lvl->rowB2, row12);
      
      /* list22(row12,row11) for each row11 in list11; keep the non-empty ones.
      ** The lists of one table are disjoint, so together they hold at most
      ** 2^width entries and laDecode22 (DECODE_CAP entries) always has room. */
      int numLists = 0;
      long used = 0;
      for (int k = 0; k < len11; ++k) {
         int n;
         const uint16_t *l22 = TABLE_ROWS(factored, table22, list11[k], &n, laDecode22 + used);
         if (!n) continue;
         lists22[numLists] = l22;
         lens22[numLists] = n; used += n; ++numLists;
      }
      work += len11;
      if (!numLists) continue;
      if (lookaheadDepth == 1) return 1;     /* --lookahead-depth 1: some row11, row12, row22 exist */
      
      for (int i13 = 0; i13 < lvl->len13; ++i13) {
         const int row13 = lvl->list13[i13];
         TableRef table23 = lvl->tables23[i13];
         if (!table23) table23 = lvl->tables23[i13] = TABLE_REF(factored, lvl->rowC2, row13);
         if (!TABLE_HAS_ROW(factored, table23, row12)) continue;
         int len23;
         const uint16_t *list23 = TABLE_ROWS(factored, table23, row12, &len23, laDecode23);
         if (lookaheadDepth == 2 && len23) return 1;   /* --lookahead-depth 2: row13, row23 also exist */
         for (int k = 0; k < len23; ++k) {
            TableRef table33 = TABLE_REF(factored, row13, list23[k]);
            for (int s = 0; s < numLists; ++s)
               for (int e = 0; e < lens22[s]; ++e)
                  if (TABLE_HAS_ROW(factored, table33, lists22[s][e])) return 1;   /* some row33 exists */
            work += used;
            if (budget >= 0 && work > budget) return -1;   /* too expensive: let the split lookahead answer */
         }
      }
   }
   return 0;
}


/* ========================================================================== *
** Split lookahead (--lookahead-method split)
**
** The lookahead asks whether a "cone" of six unknown rows exists:
**      row11, row12, row13, row22, row23 and row33.
** Every rule evaluation evolve(r1,r2,r3) = r4 splits exactly into a low part (result
** bits 0 .. lowResultBits-1, using bits 0 .. lowHalfBits-1 of r1, r2, r3) and a high
** part (result bits lowResultBits .. width-1, using bits highHalfShift .. width-1),
** which share only the two "overlap" bits highHalfShift, highHalfShift+1 of r1, r2
** and r3 (the factored half-tables store exactly these two parts).  Therefore:
**
**   a full cone exists  <=>  for some choice of the overlap bits of all six rows,
**                            a low-half cone AND a high-half cone exist with those
**                            overlap bits.
**
** The low-half cone only involves half-rows of lowHalfBits bits.  Its lists come from
** the low half-table and hold about 2 entries on average at EVERY width (the free
** overlap bits make the half problem loose), so its cost does not grow with width;
** likewise for the high half.  The full-width cone enumeration it replaces costs
** about 1.65x more per extra width (lists grow with width and are multiplied).
**
** Bookkeeping: for a given low half u11L of row11, the set of overlap vectors
** (row12, row13, row22, row23, row33; 2 bits each = 10 bits = a 1024-bit set) for
** which a low-half cone exists depends only on half-rows (see splitBuild), so it is
** kept in a per-thread cache shared by all levels (splitMemo); same for the high half.
** A candidate passes iff for some overlap value v of row11, the union of the low
** sets over the low-half row11 choices with overlap v intersects the union of the
** high sets over the high-half row11 choices with overlap v.
** -DLA_VERIFY checks every result against lookAheadReference().
** ========================================================================== */
#define OVSET_WORDS 16                 /* 1024 overlap vectors = 16 x 64 bits */

typedef struct {
   int built;                          /* data below valid for the current parent */
   int lowA, highA;                    /* half-table pair index of (rowA1,rowA2)   */
   uint64_t ctxL, ctxH;                /* half-row context keys (everything but row11) */
   int b2L, c2L, b2H, c2H;             /* halves of rowB2 / rowC2                  */
   uint8_t list12L[256], list12H[256], list13L[256], list13H[256];
   int len12L, len12H, len13L, len13H;
} SplitLevel;

static _Thread_local SplitLevel **splitLevels;
static _Thread_local int splitLevelCap;

static inline SplitLevel *splitGetLevel(int level) {
   if (__builtin_expect(level >= splitLevelCap, 0)) {
      int ncap = splitLevelCap ? splitLevelCap : 64;
      while (ncap <= level) ncap *= 2;
      splitLevels = (SplitLevel**)realloc(splitLevels, ncap * sizeof(*splitLevels));
      for (int i = splitLevelCap; i < ncap; ++i) splitLevels[i] = (SplitLevel*)calloc(1, sizeof(SplitLevel));
      splitLevelCap = ncap;
   }
   return splitLevels[level];
}

/* copy all rows of one half-table bucket (every overlap group) into out[] */
static inline int halfBucket(const uint8_t *rowsBase, const uint16_t *offs, int result, uint8_t *out) {
   const int start = offs[result << 2], end = offs[(result << 2) + 4];
   for (int i = start; i < end; i++) out[i - start] = rowsBase[i];
   return end - start;
}
#define LOW_BUCKET(pair, result, out) \
   halfBucket(lowHalfRows + ((size_t)(pair) << lowHalfBits), lowHalfOffsets + (size_t)(pair) * LOW_OFFSETS_PER_PAIR, (result), (out))
#define HIGH_BUCKET(pair, result, out) \
   halfBucket(highHalfRows + ((size_t)(pair) << highHalfBits), highHalfOffsets + (size_t)(pair) * HIGH_OFFSETS_PER_PAIR, (result), (out))

static void splitBuild(SplitLevel *sl, const row *pathRows, int candRow, int rowPhase) {
   (void)rowPhase;
   const int f = FWDOFF(rowPhase), d = DOUBLEOFF(rowPhase), t = TRIPLEOFF(rowPhase);
   const int known13 = (t >= PERIOD);
   const int a1 = pathRows[candRow - PERIOD - f], a2 = pathRows[candRow - f];
   const int b1 = pathRows[candRow - PERIOD - d], b2 = pathRows[candRow - d];
   const int c2 = pathRows[candRow - t];
   const int c1 = known13 ? pathRows[candRow + PERIOD - t] : pathRows[candRow - PERIOD - t];
   const int lb = lowHalfBits, hb = highHalfBits, hs = highHalfShift, lr = lowResultBits;
   const int lm = (1 << lb) - 1, rm = (1 << lr) - 1;
   sl->lowA  = ((a1 & lm) << lb) | (a2 & lm);
   sl->highA = ((a1 >> hs) << hb) | (a2 >> hs);
   sl->b2L = b2 & lm; sl->c2L = c2 & lm; sl->b2H = b2 >> hs; sl->c2H = c2 >> hs;
   /* The half cone depends only on these half-rows (plus row11's half): known13,
   ** rowB1, rowB2, rowA2's result bits, rowC1, rowC2.  At most 6*8+7+1 = 56 bits.   */
   {
      const int hr = highResultBits;
      uint64_t k = (uint64_t)known13;
      k = (k << lb) | (b1 & lm);  k = (k << lb) | (b2 & lm);  k = (k << lr) | (a2 & rm);
      k = (k << lb) | (c1 & lm);  k = (k << lb) | (c2 & lm);
      sl->ctxL = k;
      k = (uint64_t)known13;
      k = (k << hb) | (b1 >> hs); k = (k << hb) | (b2 >> hs); k = (k << hr) | (a2 >> lr);
      k = (k << hb) | (c1 >> hs); k = (k << hb) | (c2 >> hs);
      sl->ctxH = k;
   }
   /* row12: evolve(rowB1, rowB2, row12) = rowA2 */
   sl->len12L = LOW_BUCKET(((b1 & lm) << lb) | (b2 & lm), a2 & rm, sl->list12L);
   sl->len12H = HIGH_BUCKET(((b1 >> hs) << hb) | (b2 >> hs), a2 >> lr, sl->list12H);
   /* row13: evolve(rowC1, rowC2, row13) = rowB2, or the known row rowC1 */
   if (known13) {
      sl->list13L[0] = (uint8_t)(c1 & lm); sl->len13L = 1;
      sl->list13H[0] = (uint8_t)(c1 >> hs); sl->len13H = 1;
   } else {
      sl->len13L = LOW_BUCKET(((c1 & lm) << lb) | (c2 & lm), b2 & rm, sl->list13L);
      sl->len13H = HIGH_BUCKET(((c1 >> hs) << hb) | (c2 >> hs), b2 >> lr, sl->list13H);
   }
   sl->built = 1;
}

/* Overlap-vector set of the half cone for a fixed half of row11.  `high` selects
** the half: low halves have their overlap in the top two bits (x >> highHalfShift)
** and their result bits in the low lowResultBits bits; high halves have the overlap
** in the bottom two bits (x & 3) and their result bits above it (x >> 1).
** Set index: ov12<<8 | ov13<<6 | ov22<<4 | ov23<<2 | ov33.
** Compiled once per half (high is a compile-time constant after inlining).         */
static inline __attribute__((always_inline))
void splitHalfSetT(const SplitLevel *sl, int u11, uint64_t *set, const int high) {
   const int hs = highHalfShift, bits = high ? highHalfBits : lowHalfBits;
   const int rm = (1 << lowResultBits) - 1;
   const int resBits = high ? highResultBits : lowResultBits;
   const uint8_t *masks = high ? highHalfMasks : lowHalfMasks;
   const uint8_t *rowsAll = high ? highHalfRows : lowHalfRows;
   const uint16_t *offsAll = high ? highHalfOffsets : lowHalfOffsets;
   const size_t offsPerPair = high ? HIGH_OFFSETS_PER_PAIR : LOW_OFFSETS_PER_PAIR;
   const uint8_t *list12 = high ? sl->list12H : sl->list12L, *list13 = high ? sl->list13H : sl->list13L;
   const int len12 = high ? sl->len12H : sl->len12L, len13 = high ? sl->len13H : sl->len13L;
   const int b2 = high ? sl->b2H : sl->b2L, c2 = high ? sl->c2H : sl->c2L;
   #define OV(x)  (high ? ((x) & 3) : ((x) >> hs))
   #define RES(x) (high ? ((x) >> 1) : ((x) & rm))
   /* all rows of half-table bucket (pair, result): rows[start .. end-1] */
   #define BUCKET(pair, result, rowsOut, lenOut) do { \
      const uint16_t *o_ = offsAll + (size_t)(pair) * offsPerPair + ((result) << 2); \
      rowsOut = rowsAll + ((size_t)(pair) << bits) + o_[0]; lenOut = o_[4] - o_[0]; } while (0)
   uint8_t res22[256], sh22[256];
   for (int w = 0; w < OVSET_WORDS; w++) set[w] = 0;
   const int r11 = RES(u11);
   for (int i12 = 0; i12 < len12; i12++) {
      const int u12 = list12[i12];
      /* row22: evolve(rowB2, row12, row22) = row11 */
      const uint8_t *l22; int len22;
      BUCKET((b2 << bits) | u12, r11, l22, len22);
      if (!len22) continue;
      for (int i = 0; i < len22; i++) { res22[i] = (uint8_t)RES(l22[i]); sh22[i] = (uint8_t)(OV(l22[i]) << 4); }
      const int r12 = RES(u12), ov12 = OV(u12) << 2;
      for (int i13 = 0; i13 < len13; i13++) {
         const int u13 = list13[i13];
         /* row23: evolve(rowC2, row13, row23) = row12 */
         const uint8_t *l23; int len23;
         BUCKET((c2 << bits) | u13, r12, l23, len23);
         if (!len23) continue;
         uint64_t acc = 0;
         const uint8_t *maskRow = masks + ((size_t)(u13 << bits) << resBits);
         for (int i23 = 0; i23 < len23; i23++) {
            const int u23 = l23[i23];
            /* row33: evolve(row13, row23, row33) = row22; the mask lists its overlaps */
            const uint8_t *mask33 = maskRow + ((size_t)u23 << resBits);
            const int s23 = OV(u23) << 2;
            for (int i22 = 0; i22 < len22; i22++)
               acc |= (uint64_t)mask33[res22[i22]] << (sh22[i22] | s23);
         }
         set[ov12 | OV(u13)] |= acc;
      }
   }
   #undef OV
   #undef RES
   #undef BUCKET
}
static void splitHalfSet(const SplitLevel *sl, int high, int u11, uint64_t *set) {
   if (high) splitHalfSetT(sl, u11, set, 1);
   else      splitHalfSetT(sl, u11, set, 0);
}

/* Per-thread direct-mapped cache of half-cone sets, keyed by (context, row11 half).
** Shared by every level and node: half-rows repeat far more often than full rows.  */
#define SPLIT_CACHE_BITS 16
typedef struct { uint64_t key; uint64_t set[OVSET_WORDS]; } SplitCacheEntry;
static _Thread_local SplitCacheEntry *splitCache[2];
static inline const uint64_t *splitMemo(SplitLevel *sl, int high, int u11) {
   if (__builtin_expect(!splitCache[high], 0)) {
      splitCache[high] = (SplitCacheEntry*)malloc(sizeof(SplitCacheEntry) << SPLIT_CACHE_BITS);
      for (size_t i = 0; i < ((size_t)1 << SPLIT_CACHE_BITS); i++) splitCache[high][i].key = ~(uint64_t)0;
   }
   const uint64_t key = ((high ? sl->ctxH : sl->ctxL) << 8) | (uint64_t)u11;   /* <= 64 bits */
   uint64_t h = key * 0x9E3779B97F4A7C15ULL;
   SplitCacheEntry *e = &splitCache[high][h >> (64 - SPLIT_CACHE_BITS)];
   if (e->key != key) {
      splitHalfSet(sl, high, u11, e->set);
      e->key = key;
   }
   return e->set;
}

static void splitForgetLevel(int level) {
   if (level < splitLevelCap) splitLevels[level]->built = 0;
}
static void splitForgetAll(void) {
   for (int i = 0; i < splitLevelCap; ++i) splitLevels[i]->built = 0;
}

static int lookAheadSplit(row *pathRows, int candRow, int rowPhase) {
#ifndef QSIMPLE
   if (TRIPLEOFF(rowPhase) == PERIOD) return lookAheadReference(pathRows, candRow, rowPhase);
#elif 3*OFFSET == PERIOD
   return lookAheadReference(pathRows, candRow, rowPhase);
#endif
   SplitLevel *sl = splitGetLevel(candRow);
   if (__builtin_expect(!sl->built, 0)) splitBuild(sl, pathRows, candRow, rowPhase);
   const int cand = pathRows[candRow];
   const int lr = lowResultBits, rm = (1 << lr) - 1;
   const int candLow = cand & rm, candHigh = cand >> lr;
   /* row11: evolve(rowA1, rowA2, row11) = candidate; overlaps that both halves allow */
   int common = lowHalfMasks[((size_t)sl->lowA << lr) | candLow] & highHalfMasks[((size_t)sl->highA << highResultBits) | candHigh];
   if (!common) return 0;
   const uint16_t *offL = lowHalfOffsets + (size_t)sl->lowA * LOW_OFFSETS_PER_PAIR + (candLow << 2);
   const uint16_t *offH = highHalfOffsets + (size_t)sl->highA * HIGH_OFFSETS_PER_PAIR + (candHigh << 2);
   const uint8_t *rowsL = lowHalfRows + ((size_t)sl->lowA << lowHalfBits);
   const uint8_t *rowsH = highHalfRows + ((size_t)sl->highA << highHalfBits);
   for (int v = 0; v < 4; v++) {
      if (!((common >> v) & 1)) continue;
      uint64_t unionL[OVSET_WORDS] = {0}, unionH[OVSET_WORDS] = {0};
      for (int i = offL[v]; i < offL[v + 1]; i++) {
         const uint64_t *s = splitMemo(sl, 0, rowsL[i]);
         for (int w = 0; w < OVSET_WORDS; w++) unionL[w] |= s[w];
      }
      for (int i = offH[v]; i < offH[v + 1]; i++) {
         const uint64_t *s = splitMemo(sl, 1, rowsH[i]);
         for (int w = 0; w < OVSET_WORDS; w++) unionH[w] |= s[w];
      }
      uint64_t any = 0;
      for (int w = 0; w < OVSET_WORDS; w++) any |= unionL[w] & unionH[w];
      if (any) return 1;
   }
   return 0;
}

static int lookAheadFastPair(row *pathRows, int candRow, int rowPhase) { return lookAheadFastT(pathRows, candRow, rowPhase, 0); }
static int lookAheadFastFactored(row *pathRows, int candRow, int rowPhase) { return lookAheadFastT(pathRows, candRow, rowPhase, 1); }
static inline int lookAheadFast(row *pathRows, int candRow, int rowPhase) {
   return useFactored ? lookAheadFastFactored(pathRows, candRow, rowPhase) : lookAheadFastPair(pathRows, candRow, rowPhase);
}

/* full-depth lookahead with the method chosen for this search (cone, split or hybrid) */
static inline int lookAheadChosen(row *pathRows, int candRow, int rowPhase) {
   if (useSplitLookahead == 1) return lookAheadSplit(pathRows, candRow, rowPhase);
   if (useSplitLookahead == 2) {
      laBudget = lookaheadBudget;
      int r = lookAheadFast(pathRows, candRow, rowPhase);
      laBudget = -1;
      if (r >= 0) return r;
      return lookAheadSplit(pathRows, candRow, rowPhase);
   }
   return lookAheadFast(pathRows, candRow, rowPhase);
}

int lookAhead(row *pathRows, int candRow, int rowPhase){
#ifndef NOCACHE
   /* with the lookahead cache enabled, keep the original code path */
   #ifndef QSIMPLE
   if (params[P_CACHEMEM] != 0)
   #endif
      return lookAheadReference(pathRows, candRow, ROWPHASE_ARG);
#endif
#ifdef LA_VERIFY
   if (lookaheadDepth == 3) {
      int r1 = lookAheadChosen(pathRows, candRow, ROWPHASE_ARG);
      int r2 = lookAheadReference(pathRows, candRow, ROWPHASE_ARG);
      if (r1 != r2) { fprintf(stderr, "LOOKAHEAD MISMATCH new=%d old=%d\n", r1, r2); abort(); }
      return r1;
   }
   return lookAheadFast(pathRows, candRow, ROWPHASE_ARG);   /* weaker depth: nothing to compare */
#else
   if (lookaheadDepth == 3) return lookAheadChosen(pathRows, candRow, ROWPHASE_ARG);
   return lookAheadFast(pathRows, candRow, ROWPHASE_ARG);
#endif
}

/* Testing for subperiodic patterns
**
** For each possible phase of the ship, equivRow[0][phase] gives the row that 
** is equivalent if the pattern is subperiodic with a specified period.
** equivRow[1] is necessary if gcd(period,offset) has two distinct prime 
** divisors, as two subperiods need to be tested (e.g., if speed == 6c/12, we
** must test subperiods 4 and 6).
*/

#ifndef QSIMPLE

int equivRow[2][MAXPERIOD];

int smallestDivisor(int b){
   int c = 2;
   while (b % c) ++c;
   return c;
}

void makeEqRows(int maxFactor, int divNum){
   int tempEquivRow[MAXPERIOD];
   int i,j;
   for (i = 0; i < period; ++i){
      tempEquivRow[i] = i;
      for (j = 0; j < maxFactor; ++j){
         tempEquivRow[i] += backOff[tempEquivRow[i] % period];
      }
      tempEquivRow[i] -= offset * maxFactor + i;
      equivRow[divNum][i] = tempEquivRow[i];
   }
   for (i = 0; i < period; ++i){     /* make equivRow[i] negative if possible */
      if (tempEquivRow[i] > 0){
         equivRow[divNum][i + tempEquivRow[i]] = -1 * tempEquivRow[i];
      }
   }
}

/* make phase tables for determining equivalent subperiodic rows */
void makeSubperiodTables(void){
   if (gcd(period,offset) > 1){
      int div1 = smallestDivisor(gcd(period,offset));
      makeEqRows(period / div1,0);
      int div2 = gcd(period,offset);
      while (div2 % div1 == 0) div2 /= div1;
      if (div2 != 1)
         makeEqRows(period / smallestDivisor(div2),1);
      else                                /* If gcd(period,offset) has only one prime divisor, just  */
         makeEqRows(period / div1,1);     /* reuse it.  We don't run the subperiod check very often. */
   }
}

int subperiodTest(node x, int divNum, row *pathRows, int nodeRow, uint32_t lastRow){
   int i,a;
   node y,z;
   
   int rowPhase;
   
   rowPhase = (peekPhase(x) + lastRow - nodeRow) % period + period;
   a = lastRow;   /* lastRow == 0 when calling from queue */
   while (equivRow[divNum][rowPhase % period] >= 0){
      rowPhase--;
      a--;
   }
   
   if (pathRows != NULL){
      int b = a + equivRow[divNum][rowPhase % period];
      
      while (b > nodeRow){
         if (pathRows[a] != pathRows[b])
            return 0;
         a -= period;
         b -= period;
      }
      z = x;
      for (i=0; i < (nodeRow - b); i++) z = PARENT(z);
      
      if (a > nodeRow){
         if (pathRows[a] != ROW(z))
            return 0;
         a -= period;
         for (i=0; i < period; i++)
            z = PARENT(z);
      }
      y = x;
      for (i=0; i < (nodeRow - a); i++)
         y = PARENT(y);
   }
   else {
      y = x;
      for (i=0; i < -1 * a; i++)
         y = PARENT(y);
      z = y;
      for (i=0; i < -1 * equivRow[divNum][rowPhase % period]; i++)
         z = PARENT(z);
   }
   
   while (z != 0){
      if (ROW(y) != ROW(z))
         return 0;
      
      for (i=0; i < period; i++){
         y = PARENT(y);
         z = PARENT(z);
      }
   }
   return 1;
}

int subperiodic(node x, row *pathRows, int nodeRow, uint32_t lastRow){
   if (!params[P_FULLPERIOD] || gcd(period,offset) == 1)
      return 0;
   
   if (subperiodTest(x,0,pathRows,nodeRow,lastRow) || subperiodTest(x,1,pathRows,nodeRow,lastRow))
      return 1;
   
   return 0;
}

#endif

/* 
** process() dequeues the node at the head of the queue and enqueues any valid
** child nodes.  Spaceships are detected by a call to terminal().
*/

/* One decode buffer per search-tree level (the list of candidate rows of a level  */
/* must survive while deeper levels are explored).                                */
static _Thread_local uint16_t **levelDecodeBufs;
static _Thread_local int levelDecodeCap;
static uint16_t *levelDecodeBuf(int level) {
   if (__builtin_expect(level >= levelDecodeCap, 0)) {
      int ncap = levelDecodeCap ? levelDecodeCap : 64;
      while (ncap <= level) ncap *= 2;
      levelDecodeBufs = (uint16_t**)realloc(levelDecodeBufs, ncap * sizeof(*levelDecodeBufs));
      memset(levelDecodeBufs + levelDecodeCap, 0, (ncap - levelDecodeCap) * sizeof(*levelDecodeBufs));
      levelDecodeCap = ncap;
   }
   if (!levelDecodeBufs[level]) levelDecodeBufs[level] = (uint16_t*)malloc(DECODE_CAP * sizeof(uint16_t));
   return levelDecodeBufs[level];
}

void process(node theNode)
{
   long long int i;
   int firstRow = 0;
   int numCandidates;
   int matchFlag = 1;
   node x = theNode;
   int rowPhase = peekPhase(x);
   const row *candidates;
   row pathRows[2*MAXPERIOD + 2];
   int currRow = 2*period + 1;
   for (i = currRow - 1; i >= 0; --i){
      pathRows[i] = ROW(x);
      x = PARENT(x);
   }
   laForgetLevel(currRow);   /* new parent node: forget cached lookahead context */
   
   ++rowPhase;
   if (rowPhase == period) rowPhase = 0;
   
   getSuccessors( pathRows[currRow - 2 * PERIOD],
                   pathRows[currRow - PERIOD],
                   pathRows[currRow - PERIOD + BACKOFF(rowPhase)],
                   &candidates,
                   &numCandidates,
                   levelDecodeBuf(currRow) );
   
   /* we just ran dequeue() which changed DeepQHead so */
   /* we need to look at the previous head location    */
   uint32_t deepIndex = extensionIndex[prevExtQHead];
   
   if (theNode == 0){
      firstRow = 1;
   }
   else if (deepIndex > 1){  /* This means we have a saved extension for this node */
      
      /* Sanity check: do the extension rows match the node rows? */
      node y = theNode;
      for (i = 0; i < 2*period; ++i){
         uint16_t startRow = extensionRows[deepIndex][1] + 1;
         if (extensionRows[deepIndex][startRow - i] != ROW(y)){
            fprintf(stderr, "Warning: non-matching rows detected at node %u in process()\n",theNode);
            matchFlag = 0;
            free(extensionRows[deepIndex]);
            extensionRows[deepIndex] = 0;
            break;
         }
         y = PARENT(y);
      }
      
      if (matchFlag){
         uint16_t deepStart = extensionRows[deepIndex][1] + 2;
         ++extensionRows[deepIndex][1];
         int dropped = 0;
         
         while (firstRow < numCandidates && candidates[firstRow] != extensionRows[deepIndex][deepStart]) ++firstRow;
         if (firstRow == numCandidates){   /* saved row is not a valid successor: drop the extension */
            fprintf(stderr, "Warning: invalid extension row at node %u in process()\n", theNode);
            free(extensionRows[deepIndex]);
            extensionRows[deepIndex] = 0;
            firstRow = 0;
            dropped = 1;
         }
         else if (!isVisited(theNode, candidates[firstRow])){
            enqueue(theNode, candidates[firstRow]);
            extensionIndex[extQTail - 1] = deepIndex;
            if (currentDepth() > longest){
               if (params[P_LONGEST]) bufferPattern(qTail-1, NULL, 0, 0, 0);
               longest = currentDepth();
            }
            if (terminal(qTail-1) && !terminal(PARENT(qTail-1))) success(qTail-1, NULL, 0, 0);
            setVisited(qTail - 1);
            if (extensionRows[deepIndex][1] > extensionRows[deepIndex][0]){
               extensionIndex[extQTail - 1] = 0;
            }
         }
         else     /* flag extension for elimination if it produces a previously seen node */
            extensionRows[deepIndex][1] = extensionRows[deepIndex][0] + 1;  /* extension will be eliminated by subsequent length check */
         
         /* eliminate extension if it gets too short */
         if (!dropped && extensionRows[deepIndex][1] > extensionRows[deepIndex][0]){
            free(extensionRows[deepIndex]);
            extensionRows[deepIndex] = 0;
         }
         if (!dropped) ++firstRow;
      }
   }
   
   /* This is no longer in the queue, so we can clear it */
   extensionIndex[prevExtQHead] = 0;
   
   for (i = firstRow; i < numCandidates; ++i){
      pathRows[currRow] = candidates[i];
      if (!isVisited(theNode, pathRows[currRow]) && lookAhead(pathRows, currRow, rowPhase)){
         enqueue(theNode, pathRows[currRow]);
         if (currentDepth() > longest){
            if (params[P_LONGEST]) bufferPattern(qTail-1, NULL, 0, 0, 0);
            longest = currentDepth();
         }
         if (terminal(qTail-1) && !terminal(PARENT(qTail-1))) success(qTail-1, NULL, 0, 0);
         setVisited(qTail - 1);
      }
   }
}

/* Rebuild the depth-first state (candidate lists and positions) along a saved
** extension ext (ext[0] = last row index, ext[1] = start index, rows from ext[2]).
** Returns 1 if the extension is already at least howDeep long (nothing to do),
** -1 if a saved row is not a valid successor (corrupt extension), 0 otherwise. */
int reloadExtension(uint16_t startRow, int rowPhase, uint16_t howDeep, row *ext, const uint16_t **candEndGen, int *candLeftGen, row *pathRowsGen){
   uint16_t currRow = startRow;
   
   /* Return value if length of extension is greater than deepening amount */
   if (ext[0] >= howDeep + ext[1]) return 1;
   
   memcpy(pathRowsGen + startRow, ext + 2 + ext[1], (ext[0] - ext[1] + 1) * (long long)sizeof(*pathRowsGen));
   
   for (currRow = startRow; currRow <= startRow + (ext[0] - ext[1]); ++currRow){
      getSuccessors( pathRowsGen[currRow - 2 * PERIOD],
                      pathRowsGen[currRow - PERIOD],
                      pathRowsGen[currRow - PERIOD + BACKOFF(rowPhase)],
                      &(candEndGen[currRow]),
                      &(candLeftGen[currRow]),
                      levelDecodeBuf(currRow) );
      
      candEndGen[currRow] += candLeftGen[currRow];
      
      while (candLeftGen[currRow] > 0 && *(candEndGen[currRow] - candLeftGen[currRow]) != pathRowsGen[currRow]){
         --candLeftGen[currRow];
      }
      if (candLeftGen[currRow] == 0) return -1;
      --candLeftGen[currRow];
      
      ++rowPhase;
      if (rowPhase == period) rowPhase = 0;
   }
   return 0;
}



/* MEASUREMENT ONLY: persistent memo of failed states across deepenings */
#define MB 22
static struct { uint64_t a, b; int r; } *mT; long long mLook, mHit, mRec;
static void mKey(const row *p, int n, int ph, uint64_t *a, uint64_t *b) {
   uint64_t x = 1469598103934665603ULL ^ (uint64_t)ph, y = 0x9E3779B97F4A7C15ULL + ph;
   for (int i = 0; i < n; i++) { x = (x ^ p[i]) * 1099511628211ULL; y = (y + p[i]) * 0xff51afd7ed558ccdULL; y ^= y >> 29; }
   *a = x; *b = y;
}
static void __attribute__((destructor)) mRep(void) { fprintf(stderr, "MEMO lookups=%lld hits=%lld records=%lld\n", mLook, mHit, mRec); }
int depthFirst(node theNode, uint16_t howDeep, const uint16_t **candEnd, int *candLeft, row *pathRows, _Atomic int *remainingItems, _Atomic int *forceExit, _Atomic int *passed){
   int rowPhase = peekPhase(theNode);
   node x = theNode;
   uint32_t startRow = 2*PERIOD + 1;
   uint32_t currRow = startRow;
   
   laForgetAll();   /* new node: forget all cached lookahead contexts */
   /* depth of theNode, for recording the longest partial found while deepening */
   long nodeDepth = 1;
   for (node y = theNode; y != 0; y = PARENT(y)) nodeDepth++;
   
   int i;
   for (i = currRow - 1; i >= 0; --i){
      pathRows[i] = ROW(x);
      x = PARENT(x);
   }
   ++rowPhase;
   if (rowPhase == period) rowPhase = 0;
   
   /* Reload state if we have a previous extension */
   uint32_t extIndex = extensionIndex[extQHead + theNode - qHead];
   if (extIndex > 1){
      row *extRows;
      #pragma omp critical(findDeepIndex)
      {
         extRows = extensionRows[extIndex];
      }
      
      /* Sanity check: do the extension rows match the node rows? */
      int matchFlag = 1;
      node y = theNode;
      for (i = 0; i < 2*PERIOD; ++i){
         if (extRows[extRows[1] + 1 - i] != ROW(y)){
            fprintf(stderr, "Warning: non-matching rows detected at node %u in depthFirst()\n",theNode);
            matchFlag = 0;
            break;
         }
         y = PARENT(y);
      }
      
      int reload = matchFlag ? reloadExtension( (uint16_t) startRow,
                                                rowPhase,
                                                howDeep,
                                                extRows,
                                                candEnd,
                                                candLeft,
                                                pathRows ) : -1;
      if (reload == 1)
         return 1;   /* Return if howDeep is less than the length of the previous extension */
      
      if (reload == 0){
         currRow = startRow   + extRows[0]
                              - extRows[1]
                              + 1;
         rowPhase = (rowPhase + currRow - startRow) % period;
      }
      else if (matchFlag)
         fprintf(stderr, "Warning: invalid extension row at node %u in depthFirst()\n",theNode);
      
      /* the extension is used up (or invalid): release its slot */
      #pragma omp critical(findDeepIndex)
      {
         free(extensionRows[extIndex]);
         extensionRows[extIndex] = 0;
      }
   }
   
   extensionIndex[extQHead + theNode - qHead] = 0;
   
   getSuccessors( pathRows[currRow - 2 * PERIOD],
                   pathRows[currRow - PERIOD],
                   pathRows[currRow - PERIOD + BACKOFF(rowPhase)],
                   &(candEnd[currRow]),
                   &(candLeft[currRow]),
                   levelDecodeBuf(currRow) );
   candEnd[currRow] += candLeft[currRow];
   
   int earlyExit = MIN(params[P_NUMTHREADS], (int) (qTail - qHead)/4);
   for (;;){
      /* Back up if there are no rows left to check at this depth */
      if (!candLeft[currRow]){
         { if (!mT) mT = calloc((size_t)1 << MB, sizeof(*mT)); uint64_t a, b; mKey(pathRows + currRow - 2*PERIOD, 2*PERIOD, rowPhase, &a, &b);
           size_t i = a & ((1 << MB) - 1); int r = startRow + howDeep + 1 - currRow; mT[i].a = a; mT[i].b = b; mT[i].r = r; mRec++; }
         --currRow;
#ifndef QSIMPLE   /* The value of rowPhase doesn't matter for QSIMPLE, so avoid calculating it in the main loop. */
         if (rowPhase == 0) rowPhase = period;
         --rowPhase;
#endif
         if (currRow < startRow) return 0;
         
         continue;
      }
      
      /* Check next row with fixed-depth look ahead and add it if it passes */
      pathRows[currRow] = *(candEnd[currRow] - candLeft[currRow]);
      --candLeft[currRow];
      if (!lookAhead(pathRows, currRow, rowPhase))
         continue;
      
      ++currRow;
      laForgetLevel(currRow);   /* entering a new node: its lookahead context must be rebuilt */
      if (mT) { uint64_t a, b; int ph = rowPhase + 1; if (ph == period) ph = 0; mKey(pathRows + currRow - 2*PERIOD, 2*PERIOD, ph, &a, &b);
         size_t i = a & ((1 << MB) - 1); int r = startRow + howDeep + 1 - currRow; mLook++; if (mT[i].a == a && mT[i].b == b && mT[i].r <= r) mHit++; }
      /* Longest partial: deepening can reach further than the queue (and, if every
      ** node then dies, the search ends without the queue ever getting there).  */
      if (nodeDepth + (long)(currRow - startRow) > longest) {
         #pragma omp critical(printWhileDeepening)
         {
            if (nodeDepth + (long)(currRow - startRow) > longest) {
               if (params[P_LONGEST]) bufferPattern(theNode, pathRows, startRow - 1, currRow - 1, 0);
               longest = nodeDepth + (currRow - startRow);
            }
         }
      }
#ifndef QSIMPLE   /* The value of rowPhase doesn't matter for QSIMPLE, so avoid calculating it in the main loop. */
      ++rowPhase;
      if (rowPhase == period) rowPhase = 0;
#endif
      
      /* Test for early exit conditions */
      if ( atomic_load_explicit(forceExit, memory_order_relaxed)
           || (   params[P_EARLYEXIT]
               && atomic_load_explicit(remainingItems, memory_order_relaxed) < earlyExit
               && atomic_load_explicit(passed, memory_order_relaxed) ) )
         {
         extensionIndex[extQHead + theNode - qHead] = 1;   /* flag as success without saving extension rows */
         int earlyExitHowDeep = currRow - startRow - 1;
         if (earlyExitHowDeep >= params[P_MINEXTENSION])
            saveExtension(theNode, startRow, earlyExitHowDeep, pathRows);
         return 1;
      }
      
      /* Check if we reached the desired depth. If so,  
         check if the result is a complete spaceship */
      if (currRow > startRow + howDeep){
         /* Increment successful depth-first counter (used in early exit check) */
         atomic_fetch_add_explicit(passed, 1, memory_order_relaxed);

         /* Flag that an extension was found. This value will be changed by saveExtension() */
         extensionIndex[extQHead + theNode - qHead] = 1;
         
         /* Save the extension if it is long enough */
         if (howDeep >= params[P_MINEXTENSION]){
            saveExtension(theNode, startRow, howDeep, pathRows);
         }
         
         /* Check if the extension represents a spaceship */
         if (params[P_PRINTDEEP] == 0) return 1;
         for (i = 1; i <= PERIOD; ++i){
            if (pathRows[currRow - i]) return 1;
         }
         currRow -= PERIOD;
         for (i = 1; i <= PERIOD; ++i){
            if (causesBirth[pathRows[currRow - i]]) return 1;
         }
         
         /* If we got here, then we found a spaceship! */
         #pragma omp critical(printWhileDeepening)
         {
            success(theNode, pathRows, startRow - 1, currRow + PERIOD - 1);
         }
         if (aborting)  /* Flag for early exit if the desired number of ships has been found */
            atomic_store_explicit(forceExit, 1, memory_order_seq_cst);
         return 1;
      }
      
      /* Get the list of successor rows based on the newly added row */
      getSuccessors(pathRows[currRow - 2 * PERIOD],
                     pathRows[currRow - PERIOD],
                     pathRows[currRow - PERIOD + BACKOFF(rowPhase)],
                     &(candEnd[currRow]), &(candLeft[currRow]), levelDecodeBuf(currRow));
      candEnd[currRow] += candLeft[currRow];
   }
}

int main(int argc, char *argv[]){
   printf("%s\n",BANNER);
   
   setDefaultParams();
   parseOptions(argc, argv);
   searchSetup();
   
   /* (subperiod tables are built in searchSetup() via makeSubperiodTables()) */
   
   printf("Starting search\n");
   fflush(stdout);
   
   breadthFirst();
   
   finalReport();
   return 0;
}
