#include <map>
#include <set>
#include <iosfwd>
#include <string>

namespace stk
{
namespace io
{

typedef double Time;
typedef int Step;
typedef std::map<Time, Time> TimeContainer;
typedef std::map<Step, Step> StepContainer;

struct TolerancedTime {
    // A simple container to hold toleranced time values;
    Time min;
    Time max;
};

class IOScheduler
  {
  public:
    IOScheduler();

    IOScheduler(const IOScheduler &from);

    ~IOScheduler();

    void set_force_write(); //!< Force write on next call to scheduler

    bool write_now(Time time, Step step, Time termination_time);
    bool write_now(Time time);
    bool write_now(Step  step);

    Time adjust_dt(Time dt, Time time);

    bool add_interval(Step step, Step interval=1);
    bool add_explicit(Step step);

    bool add_interval(Time time, Time delta=0.0);
    bool add_explicit(Time time);

    bool set_lookahead(int lookahead);
    bool set_start_time(Time time);
    bool set_termination_time(Time time);
    void set_synchronize() {synchronize_ = true;}
    bool get_synchronize() {return synchronize_;}

    void set_restart_time(Time time) {restartTime_ = time;}

    bool set_signal(const std::string& signal);

    /*! Function for use by encore.  They have a use case in which
     * they reset the time of a region back to zero multiple times
     * during an execution.  Normally, the IOScheduler requires
     * that time be strictly increasing; this function resets the
     * lastTime_ value to -2.0 such that the IOScheduler does not
     * throw an error about non-increasing time.
     */
    void reset_last_time();

    void print(std::ostream &out) const;

    const std::string & name() const       {return name_;}
    void set_name(const std::string &n) {name_ = n;}

  private:

    bool force_write();

    Time next_implicit_output_time(Time time) const;
    Time next_explicit_output_time(Time time) const;

    TimeContainer::const_iterator get_time_interval(Time time, bool erase_old) const;
    TolerancedTime get_toleranced_time_range(Time time) const;

    /*! timeIntervals_ is a container of time intervals where a
     * time interval is a pair <start_time, delta_time>.  A time
     * interval specifies that starting at time 'start_time',
     * output will be written every 'delta_time' intervals.  For
     * example, the interval <0.0, 0.1> would cause output to be
     * written at 0.0, 0.1, 0.2, 0.3, ....
     *
     * The stop or end time of an interval is the begin time of
     * the following interval.  Therefore, only a single time
     * interval is active at a single time and time intervals do
     * not overlap.
     */
    mutable TimeContainer timeIntervals_;

    //! See discussion for time intervals.
    StepContainer stepIntervals_;

    std::set<Time> times_; //<! List of addtional times at which output is desired.
    std::set<Step> steps_; //<! List of addtional steps at which output is desired.

    Time tolerance_; //<! tolerance for comparing times.
    mutable Time lastTime_;  //<! Last time at which output was written.
    mutable Time firstTime_; //<! First time at which output was written.
    Time lastCalledTime_; //<! Time routine called last; to estimate dt
    mutable Step lastInterval_;

    /*!
     * Number of simulation steps to adjust timestep
     * in order to hit output time exactly. If zero, don't adjust.
     * Used by adjust_dt to calculate timestep which will hit output time.
     */
    int lookAhead_;

    mutable Time startTime_; //<! Used for backwards compatibility.

    /*! Time at which this scheduler ends.  A
     * write_now(time==terminationTime) call will return true, but
     * write_now(time >terminationTime) call will return false.
     */
    Time terminationTime_;

    /*! If this is a restarted run, the time at which it was
     *  restarted.  Used to make sure that output frequencies
     *  match between an original run and a restarted run.
     */
    Time restartTime_;

    bool forceWrite_; //<! Used to force a write
    bool synchronize_; //<! Synchronize output with outputs from other regions in this procedure.

    mutable bool initialized_; //<! True if this scheduler has been called.
    std::string name_;
  };
}
}
