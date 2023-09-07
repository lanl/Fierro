#include "profile.h"
#include <Kokkos_Core.hpp>

//
// Declaration of Event class
//
Event::Event(const std::string& name)
{
  name_ = name;
  count_ = 0;
  //total_time_ = duration_t::duration::zero();
  clock_t::time_point t1 = clock_t::now();
  total_time_ = t1 - t1;
}

void Event::start()
{
    start_time_ = clock_t::now();
}

void Event::stop()
{
    count_ += 1;
    total_time_ += (clock_t::now() - start_time_);
}

double Event::get_time_in_seconds()
{
  return std::chrono::duration_cast
         <std::chrono::nanoseconds>
         (total_time_).count() * 1.0e-9;
}

int Event::get_count()
{
    return count_;
}

std::string& Event::get_name()
{
    return name_;
}


//
// Declaration of Profile class
//
Event Profile::total = Event("total");
Event Profile::fft_forward = Event("fft_forward");
Event Profile::fft_backward = Event("fft_backward");

//
std::vector<Event*> Profile::events_;

Profile::Profile()
{
}

void Profile::start(Event & event)
{
  event.start();
}

void Profile::stop(Event & event)
{
  if (event.get_count() == 0){
    events_.push_back(&event);
  }

  event.stop();
}

void Profile::start_barrier(Event & event)
{
  Kokkos::fence();
#ifdef HAVE_CUDA
  cudaDeviceSynchronize();
#elif HAVE_OPENMP
  #pragma omp barrier
#endif
  event.start();
}

void Profile::stop_barrier(Event & event)
{
  if (event.get_count() == 0){
    events_.push_back(&event);
  }

  Kokkos::fence();
#ifdef HAVE_CUDA
  cudaDeviceSynchronize();
#elif HAVE_OPENMP
  #pragma omp barrier
#endif
  event.stop();
}

void Profile::print_one(Event &event)
{
  //printf("%s : %12.4E seconds\n", event.get_name().c_str(), event.get_time_in_seconds());
  printf("\n");
  printf("%s:\n", event.get_name().c_str());
  printf("  time: %12.4E seconds", event.get_time_in_seconds());
  printf("  count: %d", event.get_count());
  printf("  fraction: %12.4E%%", event.get_time_in_seconds()/total.get_time_in_seconds()*100.0);
  printf("\n");
}

void Profile::print()
{
  for (Event* event : events_) {
    print_one(*event);
  }
}
