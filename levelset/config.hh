#ifndef LEVELPP_CONFIG_HH
#define LEVELPP_CONFIG_HH

const int heap_init_memory=4096*32;
const int u_init_memory=65536*16;
const int temp_init_memory=65536*16;

#define GRID_CLEANING

//#define MEMORY_MESSAGES
//#define SIMPLE_SEARCH
//#define FORCE_FIRST_ORDER
//#define CORNER_BAIL_WARNING

const int max_extrapolation_fields=32;

const double tolerance=1e-12;

#endif
