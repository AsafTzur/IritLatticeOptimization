#ifndef OPTIMIZATION_CONFIG_PARSER_H
#define OPTIMIZATION_CONFIG_PARSER_H

#include "opt_config.h"

void OptimizationConfigInit(OptConfigStruct *Config);
int OptimizationConfigLoadFile(const char *FileName, OptConfigStruct *Config);
int OptimizationConfigValidate(const OptConfigStruct *Config);
void OptimizationConfigFree(OptConfigStruct *Config);

#endif
