#pragma once

#include <chrono>
#include <stdio.h>
#include <string>
#include <vector>


// Defination of Event class
class Event
{

  //using chrono_t   = std::chrono;
  using clock_t    = std::chrono::high_resolution_clock;
  using duration_t = clock_t::duration;

  private:
    clock_t::time_point  start_time_;
    duration_t           total_time_;     // total time
    int                  count_;          // current call count
    std::string          name_;           // name of event
 
  public:
    Event(const std::string& name);
    void start();
    void stop();
    double get_time_in_seconds();
    int    get_count();
    std::string& get_name();

};


// Defination of Profile class
class Profile
{

  //-----------------------------------------------
  // To use the profile class for another projects,
  // change the events. Remember to define the Events
  // in profile.cpp file since they are static members.
  public:
    static Event total;
    static Event fft_forward;
    static Event fft_backward;
  //-----------------------------------------------

  private:
    static std::vector<Event*> events_;

  public:
    Profile();
    static void start(Event &event);
    static void stop(Event &event);
    static void start_barrier(Event &event);
    static void stop_barrier(Event &event);
    static void print_one(Event &event);
    static void print();
};
