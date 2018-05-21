#undef SPMD
#undef USE_NETCDF
#undef CLM_GSWP_SUPPORT
#if ( defined ABSOFT || defined PGI )
#define PGF90
#endif
#define OFFLINE

