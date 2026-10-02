// pet_sorter_with_backscatter_omp.c
// Parallel version using OpenMP
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h>
#define TAG_SINGLES 0x82
#define TAG_TIMETAG 0x86
#define BUFFER_SIZE 4096
#define REFRESH_INTERVAL 2
#define MAX_CRYSTALS 200000
#define DEBUG 0
#define BATCH_SIZE 16  // Number of periods to process before writing to disk

static size_t debug_counter = 0;
const size_t MAX_DEBUG_EVENTS = 100;
static float (*crystal_pos)[3] = NULL;
static uint64_t total_prompt_counts[3] = {0, 0, 0};
static uint64_t total_delay_counts[3] = {0, 0, 0};
static uint64_t total_prompt_backscatter_counts[2] = {0, 0}; // [0]=type1, [1]=type2
static uint64_t total_delay_backscatter_counts[2] = {0, 0}; // [0]=type1, [1]=type2
static uint64_t total_singles_processed = 0;
static uint64_t singles_used_count = 0;

typedef struct
{
  uint16_t trID;
  uint16_t axID;
  uint8_t DOI;
  uint16_t energyBits;
  uint32_t timeBits;
  bool used;
  bool used_delay;
  uint32_t eventID;        // ← from .reference
  uint32_t sourceID;       // ← from .reference
  uint32_t comptonPhantom; // ← NEW: from .reference
} SingleEvent;

typedef struct
{
  uint32_t crystalID1;
  uint32_t crystalID2;
  float timeDiff;
  float energy1;
  float energy2;
  uint32_t type;
} CoincidenceEvent;

typedef struct
{
  float ENERGY_OFFSET, ENERGY_RESOLUTION;
  float CFD_RESOLUTION;
  float ENERGY_WINDOW_511_LOW;
  float ENERGY_WINDOW_511_HIGH;
  float ENERGY_WINDOW_307_LOW;
  float ENERGY_WINDOW_307_HIGH;
  float ENERGY_WINDOW_202_LOW;
  float ENERGY_WINDOW_202_HIGH;
  float ENERGY_WINDOW_BETA_LOW;
  float ENERGY_WINDOW_BETA_HIGH;
  float TIME_WINDOW1;
  float TIME_WINDOW2;
  float DELAY_OFFSET;
  float DELAY_WINDOW1;
  float DELAY_WINDOW2;
  float CTR_511;
  float CTR_307;
  float CTR_202;
  float sigma_CTR;
  int nbCrystalsPerRing;
  int nbCrystalsPerRsector;
  int maxCrystalRing;
  int minSectorDifference;
  char fileName_lut[256];
  char multiple_policy[32];
  size_t total_events;
  size_t processed_events;
  size_t total_cycles;
  size_t processed_cycles;
  time_t start_time;
  time_t last_update;
} Config;

// Structure to hold results for one time period
typedef struct {
  // Coincidence events
  CoincidenceEvent* events_511;
  size_t count_511;
  CoincidenceEvent* events_307;
  size_t count_307;
  CoincidenceEvent* events_202;
  size_t count_202;
  CoincidenceEvent* events_511_delay;
  size_t count_511_delay;
  CoincidenceEvent* events_307_delay;
  size_t count_307_delay;
  CoincidenceEvent* events_202_delay;
  size_t count_202_delay;
  CoincidenceEvent* events_backscatter_prompt;
  size_t count_backscatter_prompt;
  CoincidenceEvent* events_backscatter_delay;
  size_t count_backscatter_delay;
  
  // Consistency data (6 uint32_t per coincidence)
  uint32_t* consist_511;
  size_t consist_count_511;
  uint32_t* consist_307;
  size_t consist_count_307;
  uint32_t* consist_202;
  size_t consist_count_202;
  uint32_t* consist_511_delay;
  size_t consist_count_511_delay;
  uint32_t* consist_307_delay;
  size_t consist_count_307_delay;
  uint32_t* consist_202_delay;
  size_t consist_count_202_delay;
  uint32_t* consist_backscatter_prompt;
  size_t consist_count_backscatter_prompt;
  uint32_t* consist_backscatter_delay;
  size_t consist_count_backscatter_delay;
  
  // Statistics
  uint64_t prompt_counts[3];
  uint64_t delay_counts[3];
  uint64_t prompt_backscatter_counts[2];
  uint64_t delay_backscatter_counts[2];
  uint64_t singles_used;
} PeriodResult;

static Config cfg;

// Check if energy sum is within 511 keV window
static inline bool is_energy_sum_in_511_window(float energy1, float energy2)
{
    float energy_sum = energy1 + energy2;
    return (energy_sum >= cfg.ENERGY_WINDOW_511_LOW && energy_sum <= cfg.ENERGY_WINDOW_511_HIGH);
}

static inline float energy_of(const SingleEvent *e)
{
  return (float)(e->energyBits * cfg.ENERGY_RESOLUTION + cfg.ENERGY_OFFSET);
}

static inline uint32_t cid_of(const SingleEvent *e)
{
  return (uint32_t)e->axID * (uint32_t)cfg.nbCrystalsPerRing + (uint32_t)e->trID;
}

static inline bool get_crystal_xyz(uint32_t cid, float *x, float *y, float *z)
{
  if (!crystal_pos)
    return false;
  if (cid >= MAX_CRYSTALS)
    return false;
  *x = crystal_pos[cid][0];
  *y = crystal_pos[cid][1];
  *z = crystal_pos[cid][2];
  return true;
}

static inline float distance_between_crystals(uint32_t c1, uint32_t c2)
{
  float x1, y1, z1, x2, y2, z2;
  if (!get_crystal_xyz(c1, &x1, &y1, &z1))
    return -1.0f;
  if (!get_crystal_xyz(c2, &x2, &y2, &z2))
    return -1.0f;
  float dx = x1 - x2, dy = y1 - y2, dz = z1 - z2;
  return sqrtf(dx * dx + dy * dy + dz * dz);
}

static inline float expected_time_ps_from_distance(float distance_mm)
{
  const float c_mmperps = 3e-1f;
  return distance_mm / c_mmperps;
}

bool isDebugPair(SingleEvent *ev1, SingleEvent *ev2)
{
  if (ev1->trID == 601 && ev1->axID == 25 && ev2->trID == 66 && ev2->axID == 73)
    return true;
  if (ev1->trID == 784 && ev1->axID == 35 && ev2->trID == 276 && ev2->axID == 16)
    return true;
  if (ev1->trID == 746 && ev1->axID == 42 && ev2->trID == 183 && ev2->axID == 60)
    return true;
  if (ev1->trID == 635 && ev1->axID == 104 && ev2->trID == 59 && ev2->axID == 14)
    return true;
  if (ev1->trID == 893 && ev1->axID == 29 && ev2->trID == 348 && ev2->axID == 30)
    return true;
  if (ev1->trID == 170 && ev1->axID == 68 && ev2->trID == 663 && ev2->axID == 34)
    return true;
  if (ev1->trID == 684 && ev1->axID == 115 && ev2->trID == 190 && ev2->axID == 0)
    return true;
  if (ev1->trID == 449 && ev1->axID == 22 && ev2->trID == 978 && ev2->axID == 97)
    return true;
  return false;
}

void show_progress(uint64_t processed_cycles, size_t words_read, size_t total_words)
{
  double elapsed = difftime(time(NULL), cfg.start_time);
  double percent = (total_words > 0) ? (double)words_read / total_words * 100.0 : 0;
  double remaining = (elapsed > 0 && percent > 0.1) ? elapsed * (100.0 - percent) / percent : 0;
  printf("\rCycles: %" PRIu64 " | File: %.1f%% (%zu/%zu M words) | Elapsed: %02d:%02d | ETA: %02d:%02d    ",
         processed_cycles, percent,
         words_read / 1000000, total_words / 1000000,
         (int)(elapsed / 60), (int)fmod(elapsed, 60),
         (int)(remaining / 60), (int)fmod(remaining, 60));
  fflush(stdout);
}

bool isWithinFOV(uint32_t c1, uint32_t c2, float dt_sec)
{
  if (!crystal_pos)
    return false;
  if (c1 >= MAX_CRYSTALS || c2 >= MAX_CRYSTALS)
    return false;
  float x1 = crystal_pos[c1][0];
  float y1 = crystal_pos[c1][1];
  float z1 = crystal_pos[c1][2];
  float x2 = crystal_pos[c2][0];
  float y2 = crystal_pos[c2][1];
  float z2 = crystal_pos[c2][2];
  float dx = x1 - x2, dy = y1 - y2, dz = z1 - z2;
  float L_mm = sqrtf(dx * dx + dy * dy + dz * dz);
  const double c_mm_per_s = 3.0e11;
  double tau_sec = (double)L_mm / c_mm_per_s + 3.0 * (double)cfg.CTR_511;
  double dt_abs = fabs((double)dt_sec);
  return dt_abs <= tau_sec;
}

bool isValidLOR(uint32_t c1, uint32_t c2, float timeDiff, int event_type)
{
  if (!crystal_pos)
  {
    fprintf(stderr, "Crystal positions not loaded!\n");
    return false;
  }
  if (c1 >= MAX_CRYSTALS || c2 >= MAX_CRYSTALS)
  {
    fprintf(stderr, "Invalid crystal ID: %u or %u\n", c1, c2);
    return false;
  }
  float x1 = crystal_pos[c1][0];
  float y1 = crystal_pos[c1][1];
  float z1 = crystal_pos[c1][2];
  float x2 = crystal_pos[c2][0];
  float y2 = crystal_pos[c2][1];
  float z2 = crystal_pos[c2][2];
  float dx = x1 - x2;
  float dy = y1 - y2;
  float dz = z1 - z2;
  float distance_mm = sqrtf(dx * dx + dy * dy + dz * dz);
  float ctr_value;
  switch (event_type)
  {
  case 0:
    ctr_value = cfg.CTR_511;
    break;
  case 1:
    ctr_value = cfg.CTR_307;
    break;
  case 2:
    ctr_value = cfg.CTR_202;
    break;
  default:
    ctr_value = cfg.CTR_511;
  }
  float sigma_CTR = ctr_value / 2.355f;
  float delta_magin = cfg.sigma_CTR * sigma_CTR;
  float delta_magin_ps = delta_magin * 1e12;
  float c_mmperps = 3e-1;
  float expectedTime_ps = distance_mm / c_mmperps;
  float timeDiffAbs = fabsf(timeDiff);
  float timeDiffAbs_ps = timeDiffAbs * 1e12;
  bool valid = ((timeDiffAbs_ps >= expectedTime_ps - delta_magin_ps) && (timeDiffAbs_ps <= expectedTime_ps + delta_magin_ps));
  return valid;
}

bool check_sector_diff(int trID1, int trID2)
{
  int idRsector1 = trID1 / cfg.nbCrystalsPerRsector;
  int idRsector2 = trID2 / cfg.nbCrystalsPerRsector;
  int numRsectors = cfg.nbCrystalsPerRing / cfg.nbCrystalsPerRsector;
  int absDiff = abs(idRsector1 - idRsector2);
  int circularDiff = abs(absDiff - numRsectors);
  int minDiff = (absDiff < circularDiff) ? absDiff : circularDiff;
  return minDiff >= cfg.minSectorDifference;
}

bool check_ring_diff(int axID1, int axID2)
{
  int circularDiff = abs(axID1 - axID2);
  return circularDiff <= cfg.maxCrystalRing;
}

// Check for a third 511 keV event in the time window for type 2 backscatter detection
static bool has_third_511_event(SingleEvent *events, size_t n, size_t i, size_t j1, size_t j2, bool is_delay)
{
    float time_window = is_delay ? cfg.DELAY_WINDOW2 : cfg.TIME_WINDOW2;
    float time_offset = is_delay ? cfg.DELAY_OFFSET : 0.0f;
    float t1 = events[i].timeBits * cfg.CFD_RESOLUTION;
    float t2a = events[j1].timeBits * cfg.CFD_RESOLUTION;
    float t2b = events[j2].timeBits * cfg.CFD_RESOLUTION;
    
    // Determine the time window boundaries
    float window_start = fminf(t1, fminf(t2a, t2b)) - time_offset;
    float window_end = fmaxf(t1, fmaxf(t2a, t2b)) + time_window;
    
    for (size_t k = 0; k < n; k++) {
        if (k == i || k == j1 || k == j2) continue;
        
        float tk = events[k].timeBits * cfg.CFD_RESOLUTION;
        float dt = tk - window_start;
        
        if (dt < 0 || dt > (window_end - window_start)) continue;
        
        if (is_delay && events[k].used_delay) continue;
        if (!is_delay && events[k].used) continue;
        
        float ek = energy_of(&events[k]);
        if (ek >= cfg.ENERGY_WINDOW_511_LOW && ek <= cfg.ENERGY_WINDOW_511_HIGH) {
            return true;
        }
    }
    return false;
}

void parse_config(const char *path)
{
  cfg = (Config){
      .ENERGY_OFFSET = 0.05f,
      .ENERGY_RESOLUTION = 0.0025f,
      .CFD_RESOLUTION = 7e-12f,
      .ENERGY_WINDOW_511_LOW = 0.411f,
      .ENERGY_WINDOW_511_HIGH = 0.611f,
      .ENERGY_WINDOW_307_LOW = 0.275f,
      .ENERGY_WINDOW_307_HIGH = 0.355f,
      .ENERGY_WINDOW_202_LOW = 0.165f,
      .ENERGY_WINDOW_202_HIGH = 0.247f,
      .ENERGY_WINDOW_BETA_LOW = 0.150f,
      .ENERGY_WINDOW_BETA_HIGH = 0.700f,
      .TIME_WINDOW1 = 5e-9,
      .TIME_WINDOW2 = 7.5e-9,
      .DELAY_OFFSET = 500e-9,
      .DELAY_WINDOW1 = 5e-9,
      .DELAY_WINDOW2 = 7.5e-9,
      .CTR_511 = 2.5e-10,
      .CTR_307 = 2.5e-10,
      .CTR_202 = 2.5e-10,
      .sigma_CTR = 3.0f,
      .nbCrystalsPerRing = 1056,
      .nbCrystalsPerRsector = 48,
      .maxCrystalRing = 120,
      .minSectorDifference = 3,
      .fileName_lut = {0},
      .multiple_policy = "take_all",
      .total_events = 0,
      .processed_events = 0,
      .total_cycles = 0,
      .processed_cycles = 0};
  FILE *f = fopen(path, "r");
  if (!f)
  {
    perror("fail to open config file");
    exit(EXIT_FAILURE);
  }
  char line[256];
  while (fgets(line, sizeof(line), f))
  {
    float fv;
    int iv;
    char str[256];
    char *end = line + strlen(line) - 1;
    while (end > line && (*end == '\n' || *end == '\r' || *end == ';'))
    {
      *end = '\0';
      end--;
    }
    if (sscanf(line, "ENERGY_OFFSET = %f", &fv) == 1)
      cfg.ENERGY_OFFSET = fv;
    else if (sscanf(line, "ENERGY_RESOLUTION = %f", &fv) == 1)
      cfg.ENERGY_RESOLUTION = fv;
    else if (sscanf(line, "CFD_RESOLUTION = %f", &fv) == 1)
      cfg.CFD_RESOLUTION = fv;
    else if (sscanf(line, "TIME_WINDOW1 = %f", &fv) == 1)
      cfg.TIME_WINDOW1 = fv;
    else if (sscanf(line, "TIME_WINDOW2 = %f", &fv) == 1)
      cfg.TIME_WINDOW2 = fv;
    else if (sscanf(line, "DELAY_OFFSET = %f", &fv) == 1)
      cfg.DELAY_OFFSET = fv;
    else if (sscanf(line, "DELAY_WINDOW1 = %f", &fv) == 1)
      cfg.DELAY_WINDOW1 = fv;
    else if (sscanf(line, "DELAY_WINDOW2 = %f", &fv) == 1)
      cfg.DELAY_WINDOW2 = fv;
    else if (sscanf(line, "CTR_511 = %f", &fv) == 1)
      cfg.CTR_511 = fv;
    else if (sscanf(line, "CTR_307 = %f", &fv) == 1)
      cfg.CTR_307 = fv;
    else if (sscanf(line, "CTR_202 = %f", &fv) == 1)
      cfg.CTR_202 = fv;
    else if (sscanf(line, "sigma_CTR = %f", &fv) == 1)
      cfg.sigma_CTR = fv;
    else if (sscanf(line, "ENERGY_WINDOW_511_LOW = %f", &fv) == 1)
      cfg.ENERGY_WINDOW_511_LOW = fv;
    else if (sscanf(line, "ENERGY_WINDOW_511_HIGH = %f", &fv) == 1)
      cfg.ENERGY_WINDOW_511_HIGH = fv;
    else if (sscanf(line, "ENERGY_WINDOW_307_LOW = %f", &fv) == 1)
      cfg.ENERGY_WINDOW_307_LOW = fv;
    else if (sscanf(line, "ENERGY_WINDOW_307_HIGH = %f", &fv) == 1)
      cfg.ENERGY_WINDOW_307_HIGH = fv;
    else if (sscanf(line, "ENERGY_WINDOW_202_LOW = %f", &fv) == 1)
      cfg.ENERGY_WINDOW_202_LOW = fv;
    else if (sscanf(line, "ENERGY_WINDOW_202_HIGH = %f", &fv) == 1)
      cfg.ENERGY_WINDOW_202_HIGH = fv;
    else if (sscanf(line, "ENERGY_WINDOW_BETA_LOW = %f", &fv) == 1)
      cfg.ENERGY_WINDOW_BETA_LOW = fv;
    else if (sscanf(line, "ENERGY_WINDOW_BETA_HIGH = %f", &fv) == 1)
      cfg.ENERGY_WINDOW_BETA_HIGH = fv;
    else if (sscanf(line, "nbCrystalsPerRing = %d", &iv) == 1)
      cfg.nbCrystalsPerRing = iv;
    else if (sscanf(line, "nbCrystalsPerRsector = %d", &iv) == 1)
      cfg.nbCrystalsPerRsector = iv;
    else if (sscanf(line, "maxCrystalRing = %d", &iv) == 1)
      cfg.maxCrystalRing = iv;
    else if (sscanf(line, "minSectorDifference = %d", &iv) == 1)
      cfg.minSectorDifference = iv;
    else if (sscanf(line, "fileName_lut = \"%255[^\"]\"", str) == 1)
      strncpy(cfg.fileName_lut, str, sizeof(cfg.fileName_lut) - 1);
    else if (sscanf(line, "fileName_lut = \'%255[^\']\'", str) == 1)
      strncpy(cfg.fileName_lut, str, sizeof(cfg.fileName_lut) - 1);
    else if (sscanf(line, "fileName_lut = %255s", str) == 1)
      strncpy(cfg.fileName_lut, str, sizeof(cfg.fileName_lut) - 1);
    else if (sscanf(line, "multiple_policy = %31s", str) == 1)
      strncpy(cfg.multiple_policy, str, sizeof(cfg.multiple_policy) - 1);
  }
  fclose(f);
  printf("\n=== Loaded Configuration ===\n");
  printf("ENERGY_OFFSET: %.6f MeV\n", cfg.ENERGY_OFFSET);
  printf("ENERGY_RESOLUTION: %.6f MeV/bit\n", cfg.ENERGY_RESOLUTION);
  printf("CFD_RESOLUTION: %.2e s (%.1f ps)\n", cfg.CFD_RESOLUTION, cfg.CFD_RESOLUTION * 1e12);
  printf("TIME_WINDOW1: %.1e s (%.1f ns)(emission)\n", cfg.TIME_WINDOW1, cfg.TIME_WINDOW1 * 1e9);
  printf("TIME_WINDOW2: %.1e s (%.1f ns)(transmission)\n", cfg.TIME_WINDOW2, cfg.TIME_WINDOW2 * 1e9);
  printf("DELAY_OFFSET: %.1e s (%.1f ns)\n", cfg.DELAY_OFFSET, cfg.DELAY_OFFSET * 1e9);
  printf("DELAY_WINDOW1: %.1e s (%.1f ns)(emission)\n", cfg.DELAY_WINDOW1, cfg.DELAY_WINDOW1 * 1e9);
  printf("DELAY_WINDOW2: %.1e s (%.1f ns)(transmission)\n", cfg.DELAY_WINDOW2, cfg.DELAY_WINDOW2 * 1e9);
  printf("CTR_511: %.1e s (%.1f ps)\n", cfg.CTR_511, cfg.CTR_511 * 1e12);
  printf("CTR_307: %.1e s (%.1f ps)\n", cfg.CTR_307, cfg.CTR_307 * 1e12);
  printf("CTR_202: %.1e s (%.1f ps)\n", cfg.CTR_202, cfg.CTR_202 * 1e12);
  printf("sigma_CTR: %.1f\n", cfg.sigma_CTR);
  printf("511keV energy window: %.3f-%.3f MeV\n", cfg.ENERGY_WINDOW_511_LOW, cfg.ENERGY_WINDOW_511_HIGH);
  printf("307keV energy window: %.3f-%.3f MeV\n", cfg.ENERGY_WINDOW_307_LOW, cfg.ENERGY_WINDOW_307_HIGH);
  printf("202keV energy window: %.3f-%.3f MeV\n", cfg.ENERGY_WINDOW_202_LOW, cfg.ENERGY_WINDOW_202_HIGH);
  printf("nbCrystalsPerRing: %d\n", cfg.nbCrystalsPerRing);
  printf("Beta energy window: %.3f-%.3f MeV\n", cfg.ENERGY_WINDOW_BETA_LOW, cfg.ENERGY_WINDOW_BETA_HIGH);
  printf("nbCrystalsPerRsector: %d\n", cfg.nbCrystalsPerRsector);
  printf("maxCrystalRing: %d\n", cfg.maxCrystalRing);
  printf("minSectorDifference: %d\n", cfg.minSectorDifference);
  printf("multiple_policy: %s\n", cfg.multiple_policy);
  printf("LUT file: %s\n", cfg.fileName_lut[0] ? cfg.fileName_lut : "(none)");
  printf("Number of threads: %d\n", omp_get_max_threads());
  printf("============================\n");
}

void load_lut(const char *filename)
{
  FILE *f = fopen(filename, "rb");
  if (!f)
  {
    perror("Failed to open LUT file");
    exit(EXIT_FAILURE);
  }
  fseek(f, 0, SEEK_END);
  long file_size = ftell(f);
  fseek(f, 0, SEEK_SET);
  int num_crystals = file_size / (6 * sizeof(float));
  printf("Loading LUT for %d crystals...\n", num_crystals);
  crystal_pos = (float (*)[3])malloc(num_crystals * sizeof(float[3]));
  if (!crystal_pos)
  {
    perror("Memory allocation failed for crystal positions");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < num_crystals; i++)
  {
    float coords[3];
    if (fread(coords, sizeof(float), 3, f) != 3)
    {
      fprintf(stderr, "Error reading crystal %d position\n", i);
      break;
    }
    fseek(f, 3 * sizeof(float), SEEK_CUR);
    crystal_pos[i][0] = (float)coords[0];
    crystal_pos[i][1] = (float)coords[1];
    crystal_pos[i][2] = (float)coords[2];
  }
  fclose(f);
  printf("LUT loaded successfully\n");
}

bool parse_packet(unsigned char buf[8], SingleEvent *ev) {
    if (buf[0] != TAG_SINGLES) return false;

    uint64_t w = (uint64_t)buf[0] |
                 ((uint64_t)buf[1] << 8) |
                 ((uint64_t)buf[2] << 16) |
                 ((uint64_t)buf[3] << 24) |
                 ((uint64_t)buf[4] << 32) |
                 ((uint64_t)buf[5] << 40) |
                 ((uint64_t)buf[6] << 48) |
                 ((uint64_t)buf[7] << 56);
    
    ev->trID       = (w >> 8)  & 0x07FF;     // 11 bits [8:18]
    ev->axID       = (w >> 19) & 0x03FF;     // 10 bits [19:28]
    ev->DOI        = 0;                         // no DOI
    ev->energyBits = (w >> 29) & 0x01FF;     // 9 bits  [29:37]
    ev->timeBits   = (w >> 38) & 0x3FFFFFF;  // 26 bits [38:63]

    ev->used = false;
    ev->used_delay = false;
    return true;
}



//bool parse_packet(unsigned char buf[8], SingleEvent *ev) {
//    if (buf[0] != TAG_SINGLES) return false;
//    uint64_t w = (uint64_t)buf[0] | 
//	         ((uint64_t)buf[1] << 8) | 
//		 ((uint64_t)buf[2] << 16) | 
//                 ((uint64_t)buf[3] << 24) | 
//		 ((uint64_t)buf[4] << 32) | 
//		 ((uint64_t)buf[5] << 40) | 
//                 ((uint64_t)buf[6] << 48) | 
//		 ((uint64_t)buf[7] << 56);
//    ev->trID = (w >> 8) & 0x07FF;       // 11 bits [8:18]
//    ev->axID = (w >> 19) & 0x03FF;      // 10 bits [19:28]
//    ev->DOI  = (w >> 29) & 0x03;        //  2 bits [29:30]
//    ev->energyBits = (w >> 31) & 0x01FF;// 9 bits  [31:39]
//    ev->timeBits = (w >> 40) & 0x3FFFFF;// 22 bits [40:61]
//    ev->used = false;
//    ev->used_delay = false;
//    return true;
//}

static int compare_by_timeBits(const void *pa, const void *pb)
{
  const SingleEvent *a = (const SingleEvent *)pa;
  const SingleEvent *b = (const SingleEvent *)pb;
  if (a->timeBits < b->timeBits)
    return -1;
  if (a->timeBits > b->timeBits)
    return 1;
  return 0;
}

typedef struct
{
  size_t a, b;
  float Esum;
} Pair;

static size_t *collect_candidates(SingleEvent *events, size_t n, size_t i, size_t *out_count)
{
  size_t cap = 16;
  size_t count = 0;
  size_t *cand = malloc(cap * sizeof(size_t));
  if (!cand)
    return NULL;
  cand[count++] = i;
  for (size_t j = i + 1; j < n; ++j)
  {
    if (events[j].used)
      continue;
    float dt_sec = (events[j].timeBits - events[i].timeBits) * cfg.CFD_RESOLUTION;
    if (dt_sec > cfg.TIME_WINDOW1)
      break;
    float energy2 = energy_of(&events[j]);
    if (energy2 >= cfg.ENERGY_WINDOW_511_LOW && energy2 <= cfg.ENERGY_WINDOW_511_HIGH)
    {
      if (count == cap)
      {
        cap *= 2;
        size_t *tmp = realloc(cand, cap * sizeof(size_t));
        if (!tmp)
        {
          free(cand);
          return NULL;
        }
        cand = tmp;
      }
      cand[count++] = j;
    }
  }
  *out_count = count;
  return cand;
}

static size_t *collect_candidates_delay(SingleEvent *events, size_t n, size_t i, size_t *out_count)
{
  size_t cap = 16;
  size_t count = 0;
  size_t *cand = malloc(cap * sizeof(size_t));
  if (!cand)
    return NULL;
  cand[count++] = i;
  for (size_t j = i + 1; j < n; ++j)
  {
    if (events[j].used_delay)
      continue;
    float dt_sec = (events[j].timeBits - events[i].timeBits) * cfg.CFD_RESOLUTION;
    if (dt_sec < cfg.DELAY_OFFSET)
      continue;
    if (dt_sec - cfg.DELAY_OFFSET > cfg.DELAY_WINDOW1)
      break;
    float energy2 = energy_of(&events[j]);
    if (energy2 >= cfg.ENERGY_WINDOW_511_LOW && energy2 <= cfg.ENERGY_WINDOW_511_HIGH)
    {
      if (count == cap)
      {
        cap *= 2;
        size_t *tmp = realloc(cand, cap * sizeof(size_t));
        if (!tmp)
        {
          free(cand);
          return NULL;
        }
        cand = tmp;
      }
      cand[count++] = j;
    }
  }
  *out_count = count;
  return cand;
}

static Pair *form_pairs(size_t *cand, size_t cand_count, SingleEvent *events, size_t *out_pc)
{
  *out_pc = 0;
  if (cand_count < 2)
    return NULL;
  size_t max_pairs = (cand_count * (cand_count - 1)) / 2;
  Pair *pairs = malloc(max_pairs * sizeof(Pair));
  if (!pairs)
    return NULL;
  size_t pc = 0;
  for (size_t a = 0; a < cand_count; ++a)
  {
    for (size_t b = a + 1; b < cand_count; ++b)
    {
      size_t p = cand[a], q = cand[b];
      if (!check_sector_diff(events[p].trID, events[q].trID))
        continue;
      if (!check_ring_diff(events[p].axID, events[q].axID))
        continue;
      float Ep = energy_of(&events[p]);
      float Eq = energy_of(&events[q]);
      pairs[pc].a = p;
      pairs[pc].b = q;
      pairs[pc].Esum = Ep + Eq;
      pc++;
    }
  }
  if (pc == 0)
  {
    free(pairs);
    *out_pc = 0;
    return NULL;
  }
  *out_pc = pc;
  return pairs;
}

// MODIFIED: write 6 uint32_t per coincidence (including backscatter type)
static inline void write_consist_to_buffer(uint32_t **buffer, size_t *count, size_t *capacity,
                                          const SingleEvent *ev1, const SingleEvent *ev2, int backscatter_type)
{
    if (*count == *capacity) {
        *capacity *= 2;
        uint32_t *tmp = realloc(*buffer, (*capacity) * 6 * sizeof(uint32_t));
        if (!tmp) {
            fprintf(stderr, "realloc failed in write_consist_to_buffer\n");
            exit(EXIT_FAILURE);
        }
        *buffer = tmp;
    }
    
    uint32_t *ptr = *buffer + (*count) * 6;
    ptr[0] = (ev1->eventID == ev2->eventID) ? 1u : 0u;  // eid_consist
    ptr[1] = ev1->sourceID;                              // sid1
    ptr[2] = ev2->sourceID;                              // sid2
    ptr[3] = ev1->comptonPhantom;                        // cp1
    ptr[4] = ev2->comptonPhantom;                        // cp2
    ptr[5] = backscatter_type;                           // backscatter type flag
    
    (*count)++;
}

static inline void write_coincidence_to_buffer(CoincidenceEvent **buffer, size_t *count, size_t *capacity,
                                              uint32_t cid1, uint32_t cid2, float timeDiff,
                                              float energy1, float energy2, uint32_t type)
{
    if (*count == *capacity) {
        *capacity *= 2;
        CoincidenceEvent *tmp = realloc(*buffer, (*capacity) * sizeof(CoincidenceEvent));
        if (!tmp) {
            fprintf(stderr, "realloc failed in write_coincidence_to_buffer\n");
            exit(EXIT_FAILURE);
        }
        *buffer = tmp;
    }
    
    CoincidenceEvent *ptr = *buffer + *count;
    ptr->crystalID1 = cid1;
    ptr->crystalID2 = cid2;
    ptr->timeDiff = timeDiff;
    ptr->energy1 = energy1;
    ptr->energy2 = energy2;
    ptr->type = type;
    
    (*count)++;
}

static void write_511_pair_and_mark_to_buffer(SingleEvent *events, size_t p, size_t q,
                                             CoincidenceEvent **buf511, size_t *cnt511, size_t *cap511,
                                             uint32_t **consist511, size_t *ccnt511, size_t *ccap511)
{
  float dt_sec = fabsf((events[q].timeBits - events[p].timeBits) * cfg.CFD_RESOLUTION);
  uint32_t cidp = cid_of(&events[p]);
  uint32_t cidq = cid_of(&events[q]);
  
  write_coincidence_to_buffer(buf511, cnt511, cap511, cidp, cidq, (float)dt_sec, 
                             energy_of(&events[p]), energy_of(&events[q]), 0);
                             
  write_consist_to_buffer(consist511, ccnt511, ccap511, &events[p], &events[q], 0); // Not backscatter
  
  if (!events[p].used)
  {
    events[p].used = true;
  }
  if (!events[q].used)
  {
    events[q].used = true;
  }
}

static bool handle_single_pair_to_buffer(size_t p, size_t q, SingleEvent *events,
                                        CoincidenceEvent **buf511, size_t *cnt511, size_t *cap511,
                                        uint32_t **consist511, size_t *ccnt511, size_t *ccap511)
{
  float dt_sec = fabsf((events[q].timeBits - events[p].timeBits) * cfg.CFD_RESOLUTION);
  uint32_t cidp = cid_of(&events[p]);
  uint32_t cidq = cid_of(&events[q]);
  if (isWithinFOV(cidp, cidq, dt_sec))
  {
    write_511_pair_and_mark_to_buffer(events, p, q, buf511, cnt511, cap511, consist511, ccnt511, ccap511);
    return true;
  }
  return false;
}

static void handle_take_all_to_buffer(Pair *pairs, size_t pc, size_t *cand, size_t cand_count,
                                     SingleEvent *events, CoincidenceEvent **buf511, size_t *cnt511, size_t *cap511,
                                     uint32_t **consist511, size_t *ccnt511, size_t *ccap511,
                                     size_t seed_i, bool *any_written, bool *seed_used)
{
  *any_written = false;
  *seed_used = false;
  size_t max_idx = 0;
  for (size_t k = 0; k < cand_count; ++k)
    if (cand[k] > max_idx)
      max_idx = cand[k];
  bool *will_be_used = calloc(max_idx + 1, sizeof(bool));
  if (!will_be_used)
  {
    fprintf(stderr, "calloc failed\n");
    return;
  }
  for (size_t k = 0; k < pc; ++k)
  {
    size_t p = pairs[k].a, q = pairs[k].b;
    float dt_sec = fabsf((events[q].timeBits - events[p].timeBits) * cfg.CFD_RESOLUTION);
    uint32_t cidp = cid_of(&events[p]);
    uint32_t cidq = cid_of(&events[q]);
    if (isWithinFOV(cidp, cidq, dt_sec))
    {
      write_coincidence_to_buffer(buf511, cnt511, cap511, cidp, cidq, (float)dt_sec,
                                 energy_of(&events[p]), energy_of(&events[q]), 0);
      write_consist_to_buffer(consist511, ccnt511, ccap511, &events[p], &events[q], 0); // Not backscatter
      
      will_be_used[p] = true;
      will_be_used[q] = true;
      *any_written = true;
    }
  }
  if (*any_written)
  {
    for (size_t kk = 0; kk < cand_count; ++kk)
    {
      size_t ev = cand[kk];
      if (ev <= max_idx && will_be_used[ev])
      {
        if (!events[ev].used)
        {
          events[ev].used = true;
        }
        if (ev == seed_i)
          *seed_used = true;
      }
    }
  }
  free(will_be_used);
}

static void handle_take_winner_to_buffer(Pair *pairs, size_t pc, SingleEvent *events, size_t seed_i,
                                        CoincidenceEvent **buf511, size_t *cnt511, size_t *cap511,
                                        uint32_t **consist511, size_t *ccnt511, size_t *ccap511,
                                        bool *any_written, bool *seed_used)
{
  *any_written = false;
  *seed_used = false;
  bool *tried = calloc(pc, sizeof(bool));
  if (!tried)
  {
    fprintf(stderr, "calloc failed\n");
    return;
  }
  for (size_t iter = 0; iter < pc; ++iter)
  {
    size_t max_idx = SIZE_MAX;
    float max_val = -1.0f;
    for (size_t t = 0; t < pc; ++t)
    {
      if (tried[t])
        continue;
      if (pairs[t].Esum > max_val)
      {
        max_val = pairs[t].Esum;
        max_idx = t;
      }
    }
    if (max_idx == SIZE_MAX)
      break;
    size_t p = pairs[max_idx].a, q = pairs[max_idx].b;
    float dt_sec = fabsf((events[q].timeBits - events[p].timeBits) * cfg.CFD_RESOLUTION);
    uint32_t cidp = cid_of(&events[p]);
    uint32_t cidq = cid_of(&events[q]);
    if (isWithinFOV(cidp, cidq, dt_sec))
    {
      write_511_pair_and_mark_to_buffer(events, p, q, buf511, cnt511, cap511, consist511, ccnt511, ccap511);
      *any_written = true;
      if (p == seed_i || q == seed_i)
        *seed_used = true;
      break;
    }
    tried[max_idx] = true;
  }
  free(tried);
}

static void handle_kill_all(bool *any_written, bool *seed_used)
{
  *any_written = false;
  *seed_used = false;
}

static void handle_multiple_coincidences_to_buffer(size_t *cand, size_t cand_count, Pair *pairs, size_t pc,
                                                  SingleEvent *events, CoincidenceEvent **buf511, size_t *cnt511, size_t *cap511,
                                                  uint32_t **consist511, size_t *ccnt511, size_t *ccap511,
                                                  size_t seed_i, bool *any_written, bool *seed_used)
{
  if (cand_count < 3)
  {
    *any_written = false;
    *seed_used = false;
    return;
  }
  if (strncmp(cfg.multiple_policy, "kill_all", 8) == 0)
  {
    handle_kill_all(any_written, seed_used);
    return;
  }
  if (strncmp(cfg.multiple_policy, "take_all", 14) == 0)
  {
    handle_take_all_to_buffer(pairs, pc, cand, cand_count, events, buf511, cnt511, cap511, consist511, ccnt511, ccap511, seed_i, any_written, seed_used);
    return;
  }
  if (strncmp(cfg.multiple_policy, "take_winner", 10) == 0)
  {
    handle_take_winner_to_buffer(pairs, pc, events, seed_i, buf511, cnt511, cap511, consist511, ccnt511, ccap511, any_written, seed_used);
    return;
  }
  handle_take_all_to_buffer(pairs, pc, cand, cand_count, events, buf511, cnt511, cap511, consist511, ccnt511, ccap511, seed_i, any_written, seed_used);
}

static void process_511kev_coincidence_to_buffer(SingleEvent *events, size_t n, size_t i,
                                                CoincidenceEvent **buf511, size_t *cnt511, size_t *cap511,
                                                uint32_t **consist511, size_t *ccnt511, size_t *ccap511,
                                                bool *out_any_written, bool *out_seed_used)
{
  *out_any_written = false;
  *out_seed_used = false;
  float energy1 = energy_of(&events[i]);
  if (!(energy1 >= cfg.ENERGY_WINDOW_511_LOW && energy1 <= cfg.ENERGY_WINDOW_511_HIGH))
    return;
  size_t cand_count = 0;
  size_t *cand = collect_candidates(events, n, i, &cand_count);
  if (!cand)
    return;
  if (cand_count < 2)
  {
    free(cand);
    return;
  }
  size_t pc = 0;
  Pair *pairs = form_pairs(cand, cand_count, events, &pc);
  if (pc == 0)
  {
    free(cand);
    free(pairs);
    return;
  }
  if (pc == 1)
  {
    size_t p = pairs[0].a, q = pairs[0].b;
    if (handle_single_pair_to_buffer(p, q, events, buf511, cnt511, cap511, consist511, ccnt511, ccap511))
    {
      *out_any_written = true;
      if (p == i || q == i)
        *out_seed_used = true;
    }
    free(pairs);
    free(cand);
    return;
  }
  handle_multiple_coincidences_to_buffer(cand, cand_count, pairs, pc, events, buf511, cnt511, cap511, consist511, ccnt511, ccap511, i, out_any_written, out_seed_used);
  free(pairs);
  free(cand);
}

static void process_511kev_delay_coincidence_to_buffer(SingleEvent *events, size_t n, size_t i,
                                                       CoincidenceEvent **buf511_delay, size_t *cnt511_delay, size_t *cap511_delay,
                                                       uint32_t **consist511_delay, size_t *ccnt511_delay, size_t *ccap511_delay,
                                                       bool *out_any_written)
{
  if (events[i].used_delay)
  {
    return;
  }
  *out_any_written = false;
  float energy1 = energy_of(&events[i]);
  if (!(energy1 >= cfg.ENERGY_WINDOW_511_LOW && energy1 <= cfg.ENERGY_WINDOW_511_HIGH))
    return;
  size_t cand_count = 0;
  size_t *cand = collect_candidates_delay(events, n, i, &cand_count);
  if (!cand)
    return;
  if (cand_count < 2)
  {
    free(cand);
    return;
  }
  size_t seed_idx = cand[0];
  size_t partners_cap = cand_count;
  size_t partners_count = 0;
  size_t *partners = malloc(partners_cap * sizeof(size_t));
  if (!partners)
  {
    free(cand);
    return;
  }
  float *partners_timeDiff = malloc(partners_cap * sizeof(float));
  if (!partners_timeDiff)
  {
    free(cand);
    free(partners);
    return;
  }
  for (size_t jj = 1; jj < cand_count; ++jj)
  {
    size_t q = cand[jj];
    float dt_sec = (events[q].timeBits - events[seed_idx].timeBits) * cfg.CFD_RESOLUTION;
    float timeDifference = dt_sec - cfg.DELAY_OFFSET;
    if (timeDifference < 0.0f)
      continue;
    if (timeDifference > cfg.DELAY_WINDOW1)
      continue;
    if (events[q].used_delay)
      continue;
    partners[partners_count] = q;
    partners_timeDiff[partners_count] = timeDifference;
    partners_count++;
  }
  if (partners_count == 0)
  {
    free(cand);
    free(partners);
    free(partners_timeDiff);
    return;
  }
  size_t *valid_idx = malloc(partners_count * sizeof(size_t));
  float *valid_time = malloc(partners_count * sizeof(float));
  if (!valid_idx || !valid_time)
  {
    free(cand);
    free(partners);
    free(partners_timeDiff);
    free(valid_idx);
    free(valid_time);
    return;
  }
  size_t valid_count = 0;
  uint32_t cidp = cid_of(&events[seed_idx]);
  for (size_t k = 0; k < partners_count; ++k)
  {
    size_t q = partners[k];
    float timeDifference = partners_timeDiff[k];
    if (!check_sector_diff(events[seed_idx].trID, events[q].trID))
      continue;
    if (!check_ring_diff(events[seed_idx].axID, events[q].axID))
      continue;
    uint32_t cidq = cid_of(&events[q]);
    if (!isWithinFOV(cidp, cidq, timeDifference))
      continue;
    valid_idx[valid_count] = q;
    valid_time[valid_count] = timeDifference;
    valid_count++;
  }
  if (valid_count == 0)
  {
    free(cand);
    free(partners);
    free(partners_timeDiff);
    free(valid_idx);
    free(valid_time);
    return;
  }
  if (strncmp(cfg.multiple_policy, "kill_all", 8) == 0 && valid_count >= 2)
  {
    free(cand);
    free(partners);
    free(partners_timeDiff);
    free(valid_idx);
    free(valid_time);
    return;
  }
  if (strncmp(cfg.multiple_policy, "take_winner", 10) == 0 && valid_count > 1)
  {
    float maxE = -1.0f;
    size_t win = 0;
    for (size_t kk = 0; kk < valid_count; ++kk)
    {
      float E = energy_of(&events[valid_idx[kk]]);
      if (E > maxE)
      {
        maxE = E;
        win = kk;
      }
    }
    size_t q = valid_idx[win];
    float timeDifference = valid_time[win];
    uint32_t cidq = cid_of(&events[q]);
    
    write_coincidence_to_buffer(buf511_delay, cnt511_delay, cap511_delay, cidp, cidq, timeDifference,
                               energy1, energy_of(&events[q]), 0);
    write_consist_to_buffer(consist511_delay, ccnt511_delay, ccap511_delay, &events[seed_idx], &events[q], 0); // Not backscatter
    
    *out_any_written = true;
    if (!events[seed_idx].used_delay)
      events[seed_idx].used_delay = true;
    if (!events[q].used_delay)
      events[q].used_delay = true;
      
    free(cand);
    free(partners);
    free(partners_timeDiff);
    free(valid_idx);
    free(valid_time);
    return;
  }
  for (size_t kk = 0; kk < valid_count; ++kk)
  {
    size_t q = valid_idx[kk];
    float timeDifference = valid_time[kk];
    uint32_t cidq = cid_of(&events[q]);
    
    write_coincidence_to_buffer(buf511_delay, cnt511_delay, cap511_delay, cidp, cidq, timeDifference,
                               energy1, energy_of(&events[q]), 0);
    write_consist_to_buffer(consist511_delay, ccnt511_delay, ccap511_delay, &events[seed_idx], &events[q], 0); // Not backscatter
    
    *out_any_written = true;
    if (!events[seed_idx].used_delay)
      events[seed_idx].used_delay = true;
    if (!events[q].used_delay)
      events[q].used_delay = true;
  }
  free(cand);
  free(partners);
  free(partners_timeDiff);
  free(valid_idx);
  free(valid_time);
}

static void process_transmission_coincidence_mode_to_buffer(
    SingleEvent *events, size_t n, size_t i,
    bool is_delay,
    CoincidenceEvent **buf307, size_t *cnt307, size_t *cap307,
    uint32_t **consist307, size_t *ccnt307, size_t *ccap307,
    CoincidenceEvent **buf202, size_t *cnt202, size_t *cap202,
    uint32_t **consist202, size_t *ccnt202, size_t *ccap202,
    CoincidenceEvent **buf_backscatter, size_t *cnt_backscatter, size_t *cap_backscatter,
    uint32_t **consist_backscatter, size_t *ccnt_backscatter, size_t *ccap_backscatter,
    int coincidence_type)
{
  float energy1 = energy_of(&events[i]);
  uint32_t cid1 = cid_of(&events[i]);
  bool seed_in_511_window = (energy1 >= cfg.ENERGY_WINDOW_511_LOW && energy1 <= cfg.ENERGY_WINDOW_511_HIGH);
  if (is_delay)
  {
    if (seed_in_511_window && events[i].used_delay)
      return;
  }
  else
  {
    if (seed_in_511_window && events[i].used)
      return;
  }
  if (!(energy1 >= cfg.ENERGY_WINDOW_BETA_LOW && energy1 <= cfg.ENERGY_WINDOW_BETA_HIGH))
    return;
  float time_window = is_delay ? cfg.DELAY_WINDOW2 : cfg.TIME_WINDOW2;
  float time_offset = is_delay ? cfg.DELAY_OFFSET : 0.0f;
  size_t lcap307 = 16, lcap202 = 16;
  size_t p307c = 0, p202c = 0;
  size_t *p307 = malloc(lcap307 * sizeof(size_t));
  float *t307 = malloc(lcap307 * sizeof(float));
  size_t *p202 = malloc(lcap202 * sizeof(size_t));
  float *t202 = malloc(lcap202 * sizeof(float));
  if (!p307 || !t307 || !p202 || !t202)
  {
    free(p307); free(t307); free(p202); free(t202);
    return;
  }
  for (size_t j = i + 1; j < n; ++j)
  {
    if (is_delay)
    {
      if (events[j].used_delay)
        continue;
    }
    else
    {
      if (events[j].used)
        continue;
    }
    float dt_raw = (events[j].timeBits - events[i].timeBits) * cfg.CFD_RESOLUTION;
    float timeDiff = dt_raw - time_offset;
    if (timeDiff < 0.0f)
      continue;
    if (timeDiff > time_window)
      break;
    float e2 = energy_of(&events[j]);
    if (e2 >= cfg.ENERGY_WINDOW_307_LOW && e2 <= cfg.ENERGY_WINDOW_307_HIGH)
    {
      if (p307c == lcap307)
      {
        lcap307 *= 2;
        size_t *tmp_p = realloc(p307, lcap307 * sizeof(size_t));
        float  *tmp_t = realloc(t307, lcap307 * sizeof(float));
        if (!tmp_p || !tmp_t) { free(tmp_p ? tmp_p : p307); free(tmp_t ? tmp_t : t307); free(p202); free(t202); return; }
        p307 = tmp_p;
        t307 = tmp_t;
      }
      p307[p307c] = j;
      t307[p307c] = timeDiff;
      p307c++;
    }
    else if (e2 >= cfg.ENERGY_WINDOW_202_LOW && e2 <= cfg.ENERGY_WINDOW_202_HIGH)
    {
      if (p202c == lcap202)
      {
        lcap202 *= 2;
        size_t *tmp_p = realloc(p202, lcap202 * sizeof(size_t));
        float  *tmp_t = realloc(t202, lcap202 * sizeof(float));
        if (!tmp_p || !tmp_t) { free(tmp_p ? tmp_p : p202); free(tmp_t ? tmp_t : t202); free(p307); free(t307); return; }
        p202 = tmp_p;
        t202 = tmp_t;
      }
      p202[p202c] = j;
      t202[p202c] = timeDiff;
      p202c++;
    }
  }
  size_t *acc307 = malloc((p307c ? p307c : 1) * sizeof(size_t));
  size_t acc307c = 0;
  size_t *acc202 = malloc((p202c ? p202c : 1) * sizeof(size_t));
  size_t acc202c = 0;
  size_t *acc_backscatter = malloc(((p307c + p202c) ? (p307c + p202c) : 1) * sizeof(size_t));
  size_t acc_backscatterc = 0;
  if (!acc307 || !acc202 || !acc_backscatter)
  {
    free(p307); free(t307); free(p202); free(t202);
    free(acc307); free(acc202); free(acc_backscatter);
    return;
  }
  if (p307c > 0)
  {
    size_t *valid = malloc(p307c * sizeof(size_t));
    float *vtime = malloc(p307c * sizeof(float));
    size_t vcount = 0;
    for (size_t k = 0; k < p307c; ++k)
    {
      size_t q = p307[k];
      if (!check_sector_diff(events[i].trID, events[q].trID))
        continue;
      if (!check_ring_diff(events[i].axID, events[q].axID))
        continue;
      if (!isValidLOR(cid1, cid_of(&events[q]), t307[k], 1))
        continue;
      valid[vcount] = q;
      vtime[vcount] = t307[k];
      vcount++;
    }
    if (vcount > 0)
    {
      if (strncmp(cfg.multiple_policy, "kill_all", 8) == 0 && vcount >= 2)
      {
        // drop
      }
      else if (strncmp(cfg.multiple_policy, "take_winner", 10) == 0 && vcount > 1)
      {
        float maxE = -1;
        size_t win = 0;
        for (size_t k = 0; k < vcount; ++k)
        {
          float E = energy_of(&events[valid[k]]);
          if (E > maxE)
          {
            maxE = E;
            win = k;
          }
        }
        size_t q = valid[win];
        float energy2 = energy_of(&events[q]);
        // Check if this is a backscatter event (energy sum in 511 window)
        if (is_energy_sum_in_511_window(energy1, energy2)) {
            // Check for type 2 backscatter (third 511 keV event)
            int backscatter_type = has_third_511_event(events, n, i, i, q, is_delay) ? 2 : 1;
            
            write_coincidence_to_buffer(buf_backscatter, cnt_backscatter, cap_backscatter, 
                                       cid1, cid_of(&events[q]), vtime[win],
                                       energy1, energy2, 1); // type=1 for 307 keV
                                       
            write_consist_to_buffer(consist_backscatter, ccnt_backscatter, ccap_backscatter, 
                                   &events[i], &events[q], backscatter_type);
                                       
            acc_backscatter[acc_backscatterc++] = q;
        }
        else {
            // Regular 307 keV transmission
            write_coincidence_to_buffer(buf307, cnt307, cap307, 
                                       cid1, cid_of(&events[q]), vtime[win],
                                       energy1, energy2, 1);
                                       
            write_consist_to_buffer(consist307, ccnt307, ccap307, 
                                   &events[i], &events[q], 0); // Not backscatter
                                       
            acc307[acc307c++] = q;
        }
      }
      else
      {
        for (size_t k = 0; k < vcount; ++k)
        {
          size_t q = valid[k];
          float energy2 = energy_of(&events[q]);
          // Check if this is a backscatter event (energy sum in 511 window)
          if (is_energy_sum_in_511_window(energy1, energy2)) {
              // Check for type 2 backscatter (third 511 keV event)
              int backscatter_type = has_third_511_event(events, n, i, i, q, is_delay) ? 2 : 1;
              
              write_coincidence_to_buffer(buf_backscatter, cnt_backscatter, cap_backscatter, 
                                         cid1, cid_of(&events[q]), vtime[k],
                                         energy1, energy2, 1); // type=1 for 307 keV
                                         
              write_consist_to_buffer(consist_backscatter, ccnt_backscatter, ccap_backscatter, 
                                     &events[i], &events[q], backscatter_type);
                                     
              acc_backscatter[acc_backscatterc++] = q;
          }
          else {
              // Regular 307 keV transmission
              write_coincidence_to_buffer(buf307, cnt307, cap307, 
                                         cid1, cid_of(&events[q]), vtime[k],
                                         energy1, energy2, 1);
                                         
              write_consist_to_buffer(consist307, ccnt307, ccap307, 
                                     &events[i], &events[q], 0); // Not backscatter
                                     
              acc307[acc307c++] = q;
          }
        }
      }
    }
    free(valid);
    free(vtime);
  }
  if (p202c > 0)
  {
    size_t *valid = malloc(p202c * sizeof(size_t));
    float *vtime = malloc(p202c * sizeof(float));
    size_t vcount = 0;
    for (size_t k = 0; k < p202c; ++k)
    {
      size_t q = p202[k];
      if (!check_sector_diff(events[i].trID, events[q].trID))
        continue;
      if (!check_ring_diff(events[i].axID, events[q].axID))
        continue;
      if (!isValidLOR(cid1, cid_of(&events[q]), t202[k], 2))
        continue;
      valid[vcount] = q;
      vtime[vcount] = t202[k];
      vcount++;
    }
    if (vcount > 0)
    {
      if (strncmp(cfg.multiple_policy, "kill_all", 8) == 0 && vcount >= 2)
      {
        // drop
      }
      else if (strncmp(cfg.multiple_policy, "take_winner", 10) == 0 && vcount > 1)
      {
        float maxE = -1;
        size_t win = 0;
        for (size_t k = 0; k < vcount; ++k)
        {
          float E = energy_of(&events[valid[k]]);
          if (E > maxE)
          {
            maxE = E;
            win = k;
          }
        }
        size_t q = valid[win];
        float energy2 = energy_of(&events[q]);
        // Check if this is a backscatter event (energy sum in 511 window)
        if (is_energy_sum_in_511_window(energy1, energy2)) {
            // Check for type 2 backscatter (third 511 keV event)
            int backscatter_type = has_third_511_event(events, n, i, i, q, is_delay) ? 2 : 1;
            
            write_coincidence_to_buffer(buf_backscatter, cnt_backscatter, cap_backscatter, 
                                       cid1, cid_of(&events[q]), vtime[win],
                                       energy1, energy2, 2); // type=2 for 202 keV
                                       
            write_consist_to_buffer(consist_backscatter, ccnt_backscatter, ccap_backscatter, 
                                   &events[i], &events[q], backscatter_type);
                                       
            acc_backscatter[acc_backscatterc++] = q;
        }
        else {
            // Regular 202 keV transmission
            write_coincidence_to_buffer(buf202, cnt202, cap202, 
                                       cid1, cid_of(&events[q]), vtime[win],
                                       energy1, energy2, 2);
                                       
            write_consist_to_buffer(consist202, ccnt202, ccap202, 
                                   &events[i], &events[q], 0); // Not backscatter
                                   
            acc202[acc202c++] = q;
        }
      }
      else
      {
        for (size_t k = 0; k < vcount; ++k)
        {
          size_t q = valid[k];
          float energy2 = energy_of(&events[q]);
          // Check if this is a backscatter event (energy sum in 511 window)
          if (is_energy_sum_in_511_window(energy1, energy2)) {
              // Check for type 2 backscatter (third 511 keV event)
              int backscatter_type = has_third_511_event(events, n, i, i, q, is_delay) ? 2 : 1;
              
              write_coincidence_to_buffer(buf_backscatter, cnt_backscatter, cap_backscatter, 
                                         cid1, cid_of(&events[q]), vtime[k],
                                         energy1, energy2, 2); // type=2 for 202 keV
                                         
              write_consist_to_buffer(consist_backscatter, ccnt_backscatter, ccap_backscatter, 
                                     &events[i], &events[q], backscatter_type);
                                     
              acc_backscatter[acc_backscatterc++] = q;
          }
          else {
              // Regular 202 keV transmission
              write_coincidence_to_buffer(buf202, cnt202, cap202, 
                                         cid1, cid_of(&events[q]), vtime[k],
                                         energy1, energy2, 2);
                                         
              write_consist_to_buffer(consist202, ccnt202, ccap202, 
                                     &events[i], &events[q], 0); // Not backscatter
                                     
              acc202[acc202c++] = q;
          }
        }
      }
    }
    free(valid);
    free(vtime);
  }
  if (acc307c > 0 || acc202c > 0 || acc_backscatterc > 0)
  {
    if (is_delay)
    {
      if (!events[i].used_delay)
        events[i].used_delay = true;
      for (size_t k = 0; k < acc307c; ++k)
        if (!events[acc307[k]].used_delay)
          events[acc307[k]].used_delay = true;
      for (size_t k = 0; k < acc202c; ++k)
        if (!events[acc202[k]].used_delay)
          events[acc202[k]].used_delay = true;
      for (size_t k = 0; k < acc_backscatterc; ++k)
        if (!events[acc_backscatter[k]].used_delay)
          events[acc_backscatter[k]].used_delay = true;
    }
    else
    {
      if (!events[i].used)
      {
        events[i].used = true;
      }
      for (size_t k = 0; k < acc307c; ++k)
        if (!events[acc307[k]].used)
        {
          events[acc307[k]].used = true;
        }
      for (size_t k = 0; k < acc202c; ++k)
        if (!events[acc202[k]].used)
        {
          events[acc202[k]].used = true;
        }
      for (size_t k = 0; k < acc_backscatterc; ++k)
        if (!events[acc_backscatter[k]].used)
        {
          events[acc_backscatter[k]].used = true;
        }
    }
  }
  free(acc307);
  free(acc202);
  free(acc_backscatter);
  free(p307);
  free(t307);
  free(p202);
  free(t202);
}

// Process a single period and store results in PeriodResult struct
static void process_single_period(SingleEvent *events, size_t n, PeriodResult *result)
{
  // Always initialize result to zero
  memset(result, 0, sizeof(PeriodResult));

  if (n < 2) {
    free(events);
    return;
  }
  
  // Allocate initial buffers
  size_t cap511 = BUFFER_SIZE, cap307 = BUFFER_SIZE, cap202 = BUFFER_SIZE;
  size_t cap511d = BUFFER_SIZE, cap307d = BUFFER_SIZE, cap202d = BUFFER_SIZE;
  size_t cap_bs_prompt = BUFFER_SIZE, cap_bs_delay = BUFFER_SIZE;
  
  size_t ccap511 = BUFFER_SIZE, ccap307 = BUFFER_SIZE, ccap202 = BUFFER_SIZE;
  size_t ccap511d = BUFFER_SIZE, ccap307d = BUFFER_SIZE, ccap202d = BUFFER_SIZE;
  size_t ccap_bs_prompt = BUFFER_SIZE, ccap_bs_delay = BUFFER_SIZE;
  
  result->events_511 = malloc(cap511 * sizeof(CoincidenceEvent));
  result->events_307 = malloc(cap307 * sizeof(CoincidenceEvent));
  result->events_202 = malloc(cap202 * sizeof(CoincidenceEvent));
  result->events_511_delay = malloc(cap511d * sizeof(CoincidenceEvent));
  result->events_307_delay = malloc(cap307d * sizeof(CoincidenceEvent));
  result->events_202_delay = malloc(cap202d * sizeof(CoincidenceEvent));
  result->events_backscatter_prompt = malloc(cap_bs_prompt * sizeof(CoincidenceEvent));
  result->events_backscatter_delay = malloc(cap_bs_delay * sizeof(CoincidenceEvent));
  
  result->consist_511 = malloc(ccap511 * 6 * sizeof(uint32_t));
  result->consist_307 = malloc(ccap307 * 6 * sizeof(uint32_t));
  result->consist_202 = malloc(ccap202 * 6 * sizeof(uint32_t));
  result->consist_511_delay = malloc(ccap511d * 6 * sizeof(uint32_t));
  result->consist_307_delay = malloc(ccap307d * 6 * sizeof(uint32_t));
  result->consist_202_delay = malloc(ccap202d * 6 * sizeof(uint32_t));
  result->consist_backscatter_prompt = malloc(ccap_bs_prompt * 6 * sizeof(uint32_t));
  result->consist_backscatter_delay = malloc(ccap_bs_delay * 6 * sizeof(uint32_t));
  
  qsort(events, n, sizeof(SingleEvent), compare_by_timeBits);
  
  for (size_t i = 0; i < n; i++)
  {
    if (!events[i].used)
    {
      bool any_w, seed_u;
      process_511kev_coincidence_to_buffer(events, n, i, 
                                          &result->events_511, &result->count_511, &cap511,
                                          &result->consist_511, &result->consist_count_511, &ccap511,
                                          &any_w, &seed_u);
    }
    if (!events[i].used)
    {
      process_transmission_coincidence_mode_to_buffer(events, n, i, false,
                                                    &result->events_307, &result->count_307, &cap307,
                                                    &result->consist_307, &result->consist_count_307, &ccap307,
                                                    &result->events_202, &result->count_202, &cap202,
                                                    &result->consist_202, &result->consist_count_202, &ccap202,
                                                    &result->events_backscatter_prompt, &result->count_backscatter_prompt, &cap_bs_prompt,
                                                    &result->consist_backscatter_prompt, &result->consist_count_backscatter_prompt, &ccap_bs_prompt,
                                                    0);
    }
    if (!events[i].used_delay)
    {
      bool any_w;
      process_511kev_delay_coincidence_to_buffer(events, n, i,
                                                &result->events_511_delay, &result->count_511_delay, &cap511d,
                                                &result->consist_511_delay, &result->consist_count_511_delay, &ccap511d,
                                                &any_w);
    }
    if (!events[i].used_delay)
    {
      process_transmission_coincidence_mode_to_buffer(events, n, i, true,
                                                    &result->events_307_delay, &result->count_307_delay, &cap307d,
                                                    &result->consist_307_delay, &result->consist_count_307_delay, &ccap307d,
                                                    &result->events_202_delay, &result->count_202_delay, &cap202d,
                                                    &result->consist_202_delay, &result->consist_count_202_delay, &ccap202d,
                                                    &result->events_backscatter_delay, &result->count_backscatter_delay, &cap_bs_delay,
                                                    &result->consist_backscatter_delay, &result->consist_count_backscatter_delay, &ccap_bs_delay,
                                                    0);
    }
  }
  
  // Update statistics
  result->prompt_counts[0] = result->count_511;
  result->prompt_counts[1] = result->count_307;
  result->prompt_counts[2] = result->count_202;
  result->delay_counts[0] = result->count_511_delay;
  result->delay_counts[1] = result->count_307_delay;
  result->delay_counts[2] = result->count_202_delay;
  
  result->prompt_backscatter_counts[0] = 0;
  result->prompt_backscatter_counts[1] = 0;
  result->delay_backscatter_counts[0] = 0;
  result->delay_backscatter_counts[1] = 0;
  
  // Count backscatter types
  for (size_t i = 0; i < result->consist_count_backscatter_prompt; i++) {
    uint32_t bs_type = result->consist_backscatter_prompt[i * 6 + 5];
    if (bs_type == 1) result->prompt_backscatter_counts[0]++;
    else if (bs_type == 2) result->prompt_backscatter_counts[1]++;
  }
  
  for (size_t i = 0; i < result->consist_count_backscatter_delay; i++) {
    uint32_t bs_type = result->consist_backscatter_delay[i * 6 + 5];
    if (bs_type == 1) result->delay_backscatter_counts[0]++;
    else if (bs_type == 2) result->delay_backscatter_counts[1]++;
  }
  
  // Count singles used in this period (thread-safe: no global modification)
  size_t used = 0;
  for (size_t i = 0; i < n; i++) {
    if (events[i].used) used++;
  }
  result->singles_used = used;
  
  free(events);
}

// Write results from a batch of periods to output files
static void write_batch_results(PeriodResult *results, size_t batch_size,
                               FILE *fout511, FILE *fout307, FILE *fout202,
                               FILE *fout511_delay, FILE *fout307_delay, FILE *fout202_delay,
                               FILE *fout511_consist, FILE *fout307_consist, FILE *fout202_consist,
                               FILE *fout511_delay_consist, FILE *fout307_delay_consist, FILE *fout202_delay_consist,
                               FILE *fout_backscatter_prompt, FILE *fout_backscatter_prompt_consist,
                               FILE *fout_backscatter_delay, FILE *fout_backscatter_delay_consist)
{
  for (size_t i = 0; i < batch_size; i++) {
    PeriodResult *r = &results[i];
    
    // Write 511 prompt
    if (r->count_511 > 0) {
      fwrite(r->events_511, sizeof(CoincidenceEvent), r->count_511, fout511);
      fwrite(r->consist_511, sizeof(uint32_t) * 6, r->consist_count_511, fout511_consist);
    }
    
    // Write 307 prompt
    if (r->count_307 > 0) {
      fwrite(r->events_307, sizeof(CoincidenceEvent), r->count_307, fout307);
      fwrite(r->consist_307, sizeof(uint32_t) * 6, r->consist_count_307, fout307_consist);
    }
    
    // Write 202 prompt
    if (r->count_202 > 0) {
      fwrite(r->events_202, sizeof(CoincidenceEvent), r->count_202, fout202);
      fwrite(r->consist_202, sizeof(uint32_t) * 6, r->consist_count_202, fout202_consist);
    }
    
    // Write 511 delay
    if (r->count_511_delay > 0) {
      fwrite(r->events_511_delay, sizeof(CoincidenceEvent), r->count_511_delay, fout511_delay);
      fwrite(r->consist_511_delay, sizeof(uint32_t) * 6, r->consist_count_511_delay, fout511_delay_consist);
    }
    
    // Write 307 delay
    if (r->count_307_delay > 0) {
      fwrite(r->events_307_delay, sizeof(CoincidenceEvent), r->count_307_delay, fout307_delay);
      fwrite(r->consist_307_delay, sizeof(uint32_t) * 6, r->consist_count_307_delay, fout307_delay_consist);
    }
    
    // Write 202 delay
    if (r->count_202_delay > 0) {
      fwrite(r->events_202_delay, sizeof(CoincidenceEvent), r->count_202_delay, fout202_delay);
      fwrite(r->consist_202_delay, sizeof(uint32_t) * 6, r->consist_count_202_delay, fout202_delay_consist);
    }
    
    // Write backscatter prompt
    if (r->count_backscatter_prompt > 0) {
      fwrite(r->events_backscatter_prompt, sizeof(CoincidenceEvent), r->count_backscatter_prompt, fout_backscatter_prompt);
      fwrite(r->consist_backscatter_prompt, sizeof(uint32_t) * 6, r->consist_count_backscatter_prompt, fout_backscatter_prompt_consist);
    }
    
    // Write backscatter delay
    if (r->count_backscatter_delay > 0) {
      fwrite(r->events_backscatter_delay, sizeof(CoincidenceEvent), r->count_backscatter_delay, fout_backscatter_delay);
      fwrite(r->consist_backscatter_delay, sizeof(uint32_t) * 6, r->consist_count_backscatter_delay, fout_backscatter_delay_consist);
    }
    
    // Update global statistics
    {
      for (int j = 0; j < 3; j++) {
        total_prompt_counts[j] += r->prompt_counts[j];
        total_delay_counts[j] += r->delay_counts[j];
      }
      
      for (int j = 0; j < 2; j++) {
        total_prompt_backscatter_counts[j] += r->prompt_backscatter_counts[j];
        total_delay_backscatter_counts[j] += r->delay_backscatter_counts[j];
      }
      
      singles_used_count += r->singles_used;
    }
    
    // Free allocated memory
    free(r->events_511);
    free(r->events_307);
    free(r->events_202);
    free(r->events_511_delay);
    free(r->events_307_delay);
    free(r->events_202_delay);
    free(r->events_backscatter_prompt);
    free(r->events_backscatter_delay);
    
    free(r->consist_511);
    free(r->consist_307);
    free(r->consist_202);
    free(r->consist_511_delay);
    free(r->consist_307_delay);
    free(r->consist_202_delay);
    free(r->consist_backscatter_prompt);
    free(r->consist_backscatter_delay);
  }
}

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    fprintf(stderr, "usage: %s <singles.bin> <setup.config> <output_prefix>\n", argv[0]);
    return EXIT_FAILURE;
  }
  
  parse_config(argv[2]);
  cfg.start_time = time(NULL);
  cfg.last_update = cfg.start_time;
  
  if (strlen(cfg.fileName_lut) > 0)
  {
    load_lut(cfg.fileName_lut);
  }
  else
  {
    printf("Warning: No LUT file specified in config\n");
  }
  fflush(stdout);
  
  FILE *fin = fopen(argv[1], "rb");
  if (!fin)
  {
    perror("fail to open input file");
    return EXIT_FAILURE;
  }
  
  // Open .reference file
  char ref_path[256];
  char base_path[512];
  strncpy(base_path, argv[1], sizeof(base_path) - 1);
  base_path[sizeof(base_path) - 1] = '\0';
  char *dot = strrchr(base_path, '.');
  if (dot && strcmp(dot, ".dat") == 0) {
    *dot = '\0';
  }
  snprintf(ref_path, sizeof(ref_path), "%s.reference", base_path);
  FILE *fref = fopen(ref_path, "rb");
  if (!fref)
  {
    perror("Failed to open .reference file");
    return EXIT_FAILURE;
  }
  
  // Get file size for progress tracking
  fseek(fin, 0, SEEK_END);
  long file_size_bytes = ftell(fin);
  fseek(fin, 0, SEEK_SET);
  size_t total_words = (size_t)file_size_bytes / 8;
  printf("Input file: %ld bytes (%zu 8-byte words, %.2f GB)\n",
         file_size_bytes, total_words, (double)file_size_bytes / (1024.0*1024.0*1024.0));
  fflush(stdout);
  
  char out511_path[256], out307_path[256], out202_path[256];
  char out511_delay_path[256], out307_delay_path[256], out202_delay_path[256];
  char out511_consist_path[256], out307_consist_path[256], out202_consist_path[256];
  char out511_delay_consist_path[256], out307_delay_consist_path[256], out202_delay_consist_path[256];
  char out_backscatter_prompt_path[256], out_backscatter_prompt_consist_path[256];
  char out_backscatter_delay_path[256], out_backscatter_delay_consist_path[256];
  
  snprintf(out511_path, sizeof(out511_path), "%s_511_prompt.dat", argv[3]);
  snprintf(out307_path, sizeof(out307_path), "%s_307_prompt.dat", argv[3]);
  snprintf(out202_path, sizeof(out202_path), "%s_202_prompt.dat", argv[3]);
  snprintf(out511_delay_path, sizeof(out511_delay_path), "%s_511_delay.dat", argv[3]);
  snprintf(out307_delay_path, sizeof(out307_delay_path), "%s_307_delay.dat", argv[3]);
  snprintf(out202_delay_path, sizeof(out202_delay_path), "%s_202_delay.dat", argv[3]);
  snprintf(out511_consist_path, sizeof(out511_consist_path), "%s_511_prompt.consist", argv[3]);
  snprintf(out307_consist_path, sizeof(out307_consist_path), "%s_307_prompt.consist", argv[3]);
  snprintf(out202_consist_path, sizeof(out202_consist_path), "%s_202_prompt.consist", argv[3]);
  snprintf(out511_delay_consist_path, sizeof(out511_delay_consist_path), "%s_511_delay.consist", argv[3]);
  snprintf(out307_delay_consist_path, sizeof(out307_delay_consist_path), "%s_307_delay.consist", argv[3]);
  snprintf(out202_delay_consist_path, sizeof(out202_delay_consist_path), "%s_202_delay.consist", argv[3]);
  snprintf(out_backscatter_prompt_path, sizeof(out_backscatter_prompt_path), "%s_backscatter_prompt.dat", argv[3]);
  snprintf(out_backscatter_prompt_consist_path, sizeof(out_backscatter_prompt_consist_path), "%s_backscatter_prompt.consist", argv[3]);
  snprintf(out_backscatter_delay_path, sizeof(out_backscatter_delay_path), "%s_backscatter_delay.dat", argv[3]);
  snprintf(out_backscatter_delay_consist_path, sizeof(out_backscatter_delay_consist_path), "%s_backscatter_delay.consist", argv[3]);
  
  FILE *fout511 = fopen(out511_path, "wb");
  FILE *fout307 = fopen(out307_path, "wb");
  FILE *fout202 = fopen(out202_path, "wb");
  FILE *fout511_delay = fopen(out511_delay_path, "wb");
  FILE *fout307_delay = fopen(out307_delay_path, "wb");
  FILE *fout202_delay = fopen(out202_delay_path, "wb");
  FILE *fout511_consist = fopen(out511_consist_path, "wb");
  FILE *fout307_consist = fopen(out307_consist_path, "wb");
  FILE *fout202_consist = fopen(out202_consist_path, "wb");
  FILE *fout511_delay_consist = fopen(out511_delay_consist_path, "wb");
  FILE *fout307_delay_consist = fopen(out307_delay_consist_path, "wb");
  FILE *fout202_delay_consist = fopen(out202_delay_consist_path, "wb");
  FILE *fout_backscatter_prompt = fopen(out_backscatter_prompt_path, "wb");
  FILE *fout_backscatter_prompt_consist = fopen(out_backscatter_prompt_consist_path, "wb");
  FILE *fout_backscatter_delay = fopen(out_backscatter_delay_path, "wb");
  FILE *fout_backscatter_delay_consist = fopen(out_backscatter_delay_consist_path, "wb");
  
  if (!fout511 || !fout307 || !fout202 || !fout511_delay || !fout307_delay || !fout202_delay ||
      !fout511_consist || !fout307_consist || !fout202_consist ||
      !fout511_delay_consist || !fout307_delay_consist || !fout202_delay_consist ||
      !fout_backscatter_prompt || !fout_backscatter_prompt_consist ||
      !fout_backscatter_delay || !fout_backscatter_delay_consist)
  {
    perror("fail to open output files");
    return EXIT_FAILURE;
  }
  
  size_t ev_cap = 4096, cnt = 0;
  SingleEvent *events = (SingleEvent *)malloc(ev_cap * sizeof(SingleEvent));
  if (!events) { perror("malloc fail"); return EXIT_FAILURE; }
  
  uint64_t processed_cycles = 0;
  
  int num_threads = omp_get_max_threads();
  size_t batch_capacity = (size_t)num_threads * 8;
  if (batch_capacity < 64) batch_capacity = 64;
  if (batch_capacity > 16384) batch_capacity = 16384;
  
  printf("Starting processing: %d threads, batch size = %zu\n", num_threads, batch_capacity);
  fflush(stdout);
  
  typedef struct { SingleEvent *evts; size_t count; } PeriodData;
  PeriodData   *period_data    = calloc(batch_capacity, sizeof(PeriodData));
  PeriodResult *period_results = calloc(batch_capacity, sizeof(PeriodResult));
  if (!period_data || !period_results) { perror("batch malloc fail"); return EXIT_FAILURE; }
  
  size_t batch_count = 0;
  
  // ================================================================
  // BLOCK-READ MAIN LOOP
  // Old: fread(buf8,8,1,fin) = 46 BILLION fread calls for 368GB
  //      ~10ns/call overhead → ~8 minutes wasted on function calls alone
  // New: fread 8MB blocks → ~46,000 fread calls, iterate in memory
  //      Reference data also batch-read per block
  // ================================================================
  const size_t BLOCK_WORDS = 1024 * 1024;  // 1M uint64_t = 8MB
  uint64_t *block_buf = (uint64_t *)malloc(BLOCK_WORDS * sizeof(uint64_t));
  uint32_t *ref_block = (uint32_t *)malloc(BLOCK_WORDS * 3 * sizeof(uint32_t));
  if (!block_buf || !ref_block) { perror("block buf malloc fail"); return EXIT_FAILURE; }
  
  size_t total_words_read = 0;
  uint64_t total_timetags = 0;
  size_t nread;
  
  while ((nread = fread(block_buf, sizeof(uint64_t), BLOCK_WORDS, fin)) > 0)
  {
    // Pass 1: count singles in this block (fast in-memory scan)
    size_t singles_in_block = 0;
    for (size_t w = 0; w < nread; w++) {
      if ((block_buf[w] & 0xFF) == TAG_SINGLES) singles_in_block++;
    }
    
    // Batch-read all reference data for singles in this block
    if (singles_in_block > 0) {
      size_t rr = fread(ref_block, sizeof(uint32_t), singles_in_block * 3, fref);
      if (rr != singles_in_block * 3) {
        fprintf(stderr, "Error: .reference too short (need %zu, got %zu) at %" PRIu64 " singles\n",
                singles_in_block * 3, rr, total_singles_processed);
        exit(EXIT_FAILURE);
      }
    }
    
    // Pass 2: process each word
    size_t ref_idx = 0;
    for (size_t w = 0; w < nread; w++)
    {
      uint8_t tag = (uint8_t)(block_buf[w] & 0xFF);
      
      if (tag == TAG_TIMETAG)
      {
        total_timetags++;
        if (cnt >= 2)
        {
          period_data[batch_count].evts = malloc(cnt * sizeof(SingleEvent));
          if (!period_data[batch_count].evts) { perror("malloc period"); exit(EXIT_FAILURE); }
          memcpy(period_data[batch_count].evts, events, cnt * sizeof(SingleEvent));
          period_data[batch_count].count = cnt;
          batch_count++;
          
          if (batch_count == batch_capacity)
          {
            #pragma omp parallel for schedule(dynamic)
            for (size_t bi = 0; bi < batch_count; bi++) {
              process_single_period(period_data[bi].evts, period_data[bi].count, &period_results[bi]);
            }
            write_batch_results(period_results, batch_count,
                               fout511, fout307, fout202,
                               fout511_delay, fout307_delay, fout202_delay,
                               fout511_consist, fout307_consist, fout202_consist,
                               fout511_delay_consist, fout307_delay_consist, fout202_delay_consist,
                               fout_backscatter_prompt, fout_backscatter_prompt_consist,
                               fout_backscatter_delay, fout_backscatter_delay_consist);
            processed_cycles += batch_count;
            batch_count = 0;
            show_progress(processed_cycles, total_words_read + w, total_words);
          }
        }
        cnt = 0;
      }
      else if (tag == TAG_SINGLES)
      {
        uint64_t word = block_buf[w];
        SingleEvent ev;
        ev.trID       = (word >> 8)  & 0x07FF;
        ev.axID       = (word >> 19) & 0x03FF;
        
	//ev.DOI        = (word >> 29) & 0x03;
        //ev.energyBits = (word >> 31) & 0x01FF;
        //ev.timeBits   = (word >> 40) & 0x3FFFFF;
        
	ev.DOI        = 0;
	ev.energyBits = (word >> 29) & 0x01FF;
	ev.timeBits   = (word >> 38) & 0x3FFFFFF;
		
	ev.used       = false;
        ev.used_delay = false;
        ev.eventID        = ref_block[ref_idx * 3];
        ev.sourceID       = ref_block[ref_idx * 3 + 1];
        ev.comptonPhantom = ref_block[ref_idx * 3 + 2];
	
	ref_idx++;
        
        if (cnt == ev_cap) {
          ev_cap *= 2;
          SingleEvent *tmp = (SingleEvent *)realloc(events, ev_cap * sizeof(SingleEvent));
          if (!tmp) { perror("realloc fail"); exit(EXIT_FAILURE); }
          events = tmp;
        }
        events[cnt++] = ev;
        total_singles_processed++;
      }
    }
    total_words_read += nread;
  }
  
  // Leftover singles not delimited by trailing time tag
  if (cnt >= 2) {
    period_data[batch_count].evts = malloc(cnt * sizeof(SingleEvent));
    if (!period_data[batch_count].evts) { perror("malloc"); exit(EXIT_FAILURE); }
    memcpy(period_data[batch_count].evts, events, cnt * sizeof(SingleEvent));
    period_data[batch_count].count = cnt;
    batch_count++;
  }
  
  // Flush remaining batch
  if (batch_count > 0) {
    #pragma omp parallel for schedule(dynamic)
    for (size_t bi = 0; bi < batch_count; bi++) {
      process_single_period(period_data[bi].evts, period_data[bi].count, &period_results[bi]);
    }
    write_batch_results(period_results, batch_count,
                       fout511, fout307, fout202,
                       fout511_delay, fout307_delay, fout202_delay,
                       fout511_consist, fout307_consist, fout202_consist,
                       fout511_delay_consist, fout307_delay_consist, fout202_delay_consist,
                       fout_backscatter_prompt, fout_backscatter_prompt_consist,
                       fout_backscatter_delay, fout_backscatter_delay_consist);
    processed_cycles += batch_count;
  }
  
  // Cleanup
  free(events);
  free(block_buf);
  free(ref_block);
  free(period_data);
  free(period_results);
  fclose(fin);
  fclose(fref);
  fclose(fout511); fclose(fout307); fclose(fout202);
  fclose(fout511_delay); fclose(fout307_delay); fclose(fout202_delay);
  fclose(fout511_consist); fclose(fout307_consist); fclose(fout202_consist);
  fclose(fout511_delay_consist); fclose(fout307_delay_consist); fclose(fout202_delay_consist);
  fclose(fout_backscatter_prompt); fclose(fout_backscatter_delay);
  fclose(fout_backscatter_prompt_consist); fclose(fout_backscatter_delay_consist);
  if (crystal_pos) free(crystal_pos);
  
  time_t end_time = time(NULL);
  double total_sec = difftime(end_time, cfg.start_time);
  int hrs = (int)(total_sec / 3600);
  int mins = (int)(fmod(total_sec, 3600) / 60);
  int secs = (int)fmod(total_sec, 60);
  printf("\nProcessing completed in %02d:%02d:%02d\n", hrs, mins, secs);
  printf("\n===== Summary =====\n");
  printf("Total time tags in file: %" PRIu64 "\n", total_timetags);
  printf("Total cycles processed (>=2 singles): %" PRIu64 "\n", processed_cycles);
  printf("Total singles processed: %" PRIu64 "\n", total_singles_processed);
  printf("Total singles marked used: %" PRIu64 "\n", singles_used_count);
  printf("Prompt coincidences: 511=%" PRIu64 "  307=%" PRIu64 "  202=%" PRIu64 "\n",
         total_prompt_counts[0], total_prompt_counts[1], total_prompt_counts[2]);
  printf("Delay  coincidences: 511=%" PRIu64 "  307=%" PRIu64 "  202=%" PRIu64 "\n",
         total_delay_counts[0], total_delay_counts[1], total_delay_counts[2]);
  printf("Prompt backscatter: Type1=%" PRIu64 "  Type2=%" PRIu64 "\n",
         total_prompt_backscatter_counts[0], total_prompt_backscatter_counts[1]);
  printf("Delay  backscatter: Type1=%" PRIu64 "  Type2=%" PRIu64 "\n",
         total_delay_backscatter_counts[0], total_delay_backscatter_counts[1]);
  printf("====================\n");
  return EXIT_SUCCESS;
}
