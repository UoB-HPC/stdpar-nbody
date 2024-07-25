#pragma once

using clock_timer = std::chrono::steady_clock;
using dur_t       = std::chrono::duration<double>;

template <typename F>
dur_t time(F&& f) {
  auto s = clock_timer::now();
  f();
  return dur_t(clock_timer::now() - s);
}
