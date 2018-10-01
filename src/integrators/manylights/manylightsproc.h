#ifndef MANYLIGHTSPROC_H_
#define MANYLIGHTSPROC_H_

#include "mitsuba/core/sched.h"

#include <string>

class ManyLightsWorkUnit : public WorkUnit{
public:
    void set(const WorkUnit *other);
    void load(Stream *stream);
    void save(Stream *stream) const;
    std::string toString() const;

    MTS_DECLARE_CLASS()

private:

};

MTS_IMPLEMENT_CLASS(ManyLightsWorkUnit, false, WorkUnit)

class ManyLightsWorkResult : public WorkResult{
public:
    void load(Stream *stream);
    void save(Stream *stream) const;
    void std::string toString() const;

    MTS_DECLARE_CLASS();


};

MTS_IMPLEMENT_CLASS(ManyLightsWorkResult, false, WorkResult)

class ManyLightsWorkProcessor : public WorkProcessor{
public:
    ManyLightsWorkProcessor();
    ManyLightsWorkProcessor(Stream *stream, InstanceManager *manager);
    ref<WorkUnit> createWorkUnit();
    ref<WorkResult> createWorkResult();
    ref<WorkProcessor> clone();ManyLightsWorkProcessor();
    void prepare();
    void process(const WorkUnit *work_unit, WorkResult *work_result, const bool& stop);

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS(ManyLightsWorkProcessor, false, WorkProcessor)

class ManyLightsParallelProcess : public ParallelProcess{
public:
    ManyLightsParallelProcess();

    EStatus generateWork(WorkUnit *unit, int worker);
    void processResult(const WorkResult *result, bool cancelled);
    ref<WorkProcessor> createWorkProcessor() const;

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS(ManyLightsWorkResult, false, WorkResult)

#endif