#ifndef MISCELLANEOUS_HEADER
#define MISCELLANEOUS_HEADER

#include <stdio.h>
#include "global_struct.h"

  void coarsest_level_resets( level_struct* l, struct Thread* threading );
  void set_some_coarsest_level_improvs_params_for_setup( level_struct* l, struct Thread* threading );
  void set_some_coarsest_level_improvs_params_for_solve( level_struct* l, struct Thread* threading );

#endif
