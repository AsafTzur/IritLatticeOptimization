#ifndef LEGACY_FEA_H
#define LEGACY_FEA_H

struct OptWorkflowStateStruct;

/* Legacy solver adapter.  This owns the CalculiX backend and the temporary
 * hard-coded ANSYS INP analysis additions. */
static int RunSimulation(struct OptWorkflowStateStruct *Opt);
static int AddAnsysAnalysisDefinition(struct OptWorkflowStateStruct *Opt);
static int AddAnsysAnalysisDefinitionWinglet(
    struct OptWorkflowStateStruct *Opt);
static int RunAnsysBatch(struct OptWorkflowStateStruct *Opt);
static int RunCalculixBatch(struct OptWorkflowStateStruct *Opt);

#endif /* LEGACY_FEA_H */
