#ifndef LEGACY_FEA_H
#define LEGACY_FEA_H

struct OptWorkflowStateStruct;

/* Legacy solver adapter.  This owns the CalculiX backend and the temporary
 * configuration-driven ANSYS INP analysis additions. */
static int RunSimulation(struct OptWorkflowStateStruct *Opt);
static int AddAnsysStructuralDefinition(
    struct OptWorkflowStateStruct *Opt);
static int RunAnsysBatch(struct OptWorkflowStateStruct *Opt);
static int RunCalculixBatch(struct OptWorkflowStateStruct *Opt);

#endif /* LEGACY_FEA_H */
