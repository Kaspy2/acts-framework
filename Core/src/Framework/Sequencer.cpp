// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Framework/Sequencer.hpp"

#include <algorithm>
#include <chrono>
#include <exception>
#include <numeric>

#include <TROOT.h>
#include <dfe/dfe_namedtuple.hpp>
#include <tbb/tbb.h>

#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Paths.hpp"

FW::Sequencer::Sequencer(const Sequencer::Config& cfg)
  : m_cfg(cfg), m_logger(Acts::getDefaultLogger("Sequencer", m_cfg.logLevel))
{
  // automatically determine the number of concurrent threads to use
  if (m_cfg.numThreads < 0) {
    m_cfg.numThreads = tbb::task_scheduler_init::default_num_threads();
  }
  ROOT::EnableThreadSafety();
}

void
FW::Sequencer::addService(std::shared_ptr<IService> service)
{
  if (not service) {
    throw std::invalid_argument("Can not add empty/NULL service");
  }
  m_services.push_back(std::move(service));
  ACTS_INFO("Added service '" << m_services.back()->name() << "'");
}

void
FW::Sequencer::addContextDecorator(std::shared_ptr<IContextDecorator> decorator)
{
  if (not decorator) {
    throw std::invalid_argument("Can not add empty/NULL context decorator");
  }
  m_decorators.push_back(std::move(decorator));
  ACTS_INFO("Added context decarator '" << m_decorators.back()->name() << "'");
}

void
FW::Sequencer::addReader(std::shared_ptr<IReader> reader)
{
  if (not reader) {
    throw std::invalid_argument("Can not add empty/NULL reader");
  }
  m_readers.push_back(std::move(reader));
  ACTS_INFO("Added reader '" << m_readers.back()->name() << "'");
}

void
FW::Sequencer::addAlgorithm(std::shared_ptr<IAlgorithm> algorithm)
{
  if (not algorithm) {
    throw std::invalid_argument("Can not add empty/NULL algorithm");
  }
  m_algorithms.push_back(std::move(algorithm));
  ACTS_INFO("Added algorithm '" << m_algorithms.back()->name() << "'");
}

void
FW::Sequencer::addWriter(std::shared_ptr<IWriter> writer)
{
  if (not writer) {
    throw std::invalid_argument("Can not add empty/NULL writer");
  }
  m_writers.push_back(std::move(writer));
  ACTS_INFO("Added writer '" << m_writers.back()->name() << "'");
}

std::vector<std::string>
FW::Sequencer::listAlgorithmNames() const
{
  std::vector<std::string> names;

  // WARNING this must be done in the same order as in the processing
  for (const auto& service : m_services) {
    names.push_back("Service:" + service->name());
  }
  for (const auto& decorator : m_decorators) {
    names.push_back("Decorator:" + decorator->name());
  }
  for (const auto& reader : m_readers) {
    names.push_back("Reader:" + reader->name());
  }
  for (const auto& algorithm : m_algorithms) {
    names.push_back("Algorithm:" + algorithm->name());
  }
  for (const auto& writer : m_writers) {
    names.push_back("Writer:" + writer->name());
  }

  return names;
}

std::size_t
FW::Sequencer::determineEndEvent() const
{
  // beware of possible overflow due to events == SIZE_MAX
  // use saturated add from http://locklessinc.com/articles/sat_arithmetic/
  std::size_t endRequested = m_cfg.skip + m_cfg.events;
  endRequested |= -(endRequested < m_cfg.events);

  // determine maximum events available from readers
  std::size_t endOnFile      = SIZE_MAX;
  auto        shortestReader = std::min_element(
      m_readers.begin(), m_readers.end(), [](const auto& a, const auto& b) {
        return (a->numEvents() < b->numEvents());
      });
  if (shortestReader != m_readers.end()) {
    endOnFile = (*shortestReader)->numEvents();
  }

  // without readers, skipping has no meaning
  if (m_readers.empty() and (0 < m_cfg.skip)) {
    ACTS_ERROR("Can not skip events without configured readers");
    return SIZE_MAX;
  }
  // trying to skip too many events must be an error
  if (endOnFile <= m_cfg.skip) {
    ACTS_ERROR("Less events available than requested to skip");
    return SIZE_MAX;
  }
  std::size_t endEvent = std::min(endRequested, endOnFile);

  if (endEvent == SIZE_MAX) {
    ACTS_ERROR("Could not determine number of events");
    return SIZE_MAX;
  }
  if (endOnFile < endRequested) {
    ACTS_INFO("Restrict number of events to available events");
  }

  return endEvent;
}

// helpers for per-algorithm timing information
namespace {

using Clock       = std::chrono::high_resolution_clock;
using Duration    = Clock::duration;
using Timepoint   = Clock::time_point;
using Seconds     = std::chrono::duration<double>;
using NanoSeconds = std::chrono::duration<double, std::nano>;

// RAII-based stopwatch to time execution within a block
struct StopWatch
{
  Timepoint start;
  Duration& store;

  StopWatch(Duration& s) : start(Clock::now()), store(s) {}
  ~StopWatch() { store += Clock::now() - start; }
};

// Convert duration to a printable string w/ reasonable unit.
template <typename D>
inline std::string
asString(D duration)
{
  double ns = std::chrono::duration_cast<NanoSeconds>(duration).count();
  if (1e9 < std::abs(ns)) {
    return std::to_string(ns / 1e9) + " s";
  } else if (1e6 < std::abs(ns)) {
    return std::to_string(ns / 1e6) + " ms";
  } else if (1e3 < std::abs(ns)) {
    return std::to_string(ns / 1e3) + " us";
  } else {
    return std::to_string(ns) + " ns";
  }
}

// Convert duration scaled to one event to a printable string.
template <typename D>
inline std::string
perEvent(D duration, size_t numEvents)
{
  return asString(duration / numEvents) + "/event";
}

// Store timing data
struct TimingInfo
{
  std::string identifier;
  double      time_total_s;
  double      time_perevent_s;

  DFE_NAMEDTUPLE(TimingInfo, identifier, time_total_s, time_perevent_s);
};
void
storeTiming(const std::vector<std::string>& identifiers,
            const std::vector<Duration>&    durations,
            std::size_t                     numEvents,
            std::string                     path)
{
  dfe::TsvNamedTupleWriter<TimingInfo> writer(std::move(path), 4);
  for (size_t i = 0; i < identifiers.size(); ++i) {
    TimingInfo info;
    info.identifier = identifiers[i];
    info.time_total_s
        = std::chrono::duration_cast<Seconds>(durations[i]).count();
    info.time_perevent_s = info.time_total_s / numEvents;
    writer.append(info);
  }
}

}  // namespace

int
FW::Sequencer::run()
{
  // measure overall wall clock
  Timepoint clockWallStart = Clock::now();

  // processing only works w/ a well-known number of events
  // error message is already handled by the helper function
  std::size_t endEvent = determineEndEvent();
  if (endEvent == SIZE_MAX) { return EXIT_FAILURE; }

  ACTS_INFO("Starting event loop with " << m_cfg.numThreads << " threads");
  ACTS_INFO("  " << m_services.size() << " services");
  ACTS_INFO("  " << m_decorators.size() << " context decorators");
  ACTS_INFO("  " << m_readers.size() << " readers");
  ACTS_INFO("  " << m_algorithms.size() << " algorithms");
  ACTS_INFO("  " << m_writers.size() << " writers");

  std::vector<std::string> names = listAlgorithmNames();
  std::vector<Duration>    clocksAlgorithms(names.size(), Duration::zero());
  tbb::queuing_mutex       clocksAlgorithmsMutex;

  // run start-of-run hooks
  for (auto& service : m_services) {
    names.push_back("Service:" + service->name() + ":startRun");
    clocksAlgorithms.push_back(Duration::zero());
    StopWatch sw(clocksAlgorithms.back());
    service->startRun();
  }

  // execute the parallel event loop
  tbb::task_scheduler_init init(m_cfg.numThreads);
  tbb::parallel_for(
      tbb::blocked_range<size_t>(m_cfg.skip, endEvent),
      [&](const tbb::blocked_range<size_t>& r) {
        std::vector<Duration> localClocksAlgorithms(names.size(),
                                                    Duration::zero());

        for (size_t event = r.begin(); event != r.end(); ++event) {
          // Use per-event store
          WhiteBoard eventStore(Acts::getDefaultLogger(
              "EventStore#" + std::to_string(event), m_cfg.logLevel));
          // If we ever wanted to run algorithms in parallel, this needs to be
          // changed to Algorithm context copies
          AlgorithmContext context(0, event, eventStore);
          size_t           ialgo = 0;

          // Prepare event store w/ service information
          for (auto& service : m_services) {
            StopWatch sw(localClocksAlgorithms[ialgo++]);
            service->prepare(++context);
          }
          /// Decorate the context
          for (auto& cdr : m_decorators) {
            StopWatch sw(localClocksAlgorithms[ialgo++]);
            if (cdr->decorate(++context) != ProcessCode::SUCCESS) {
              throw std::runtime_error("Failed to decorate event context");
            }
          }
          // Read everything in
          for (auto& rdr : m_readers) {
            StopWatch sw(localClocksAlgorithms[ialgo++]);
            if (rdr->read(++context) != ProcessCode::SUCCESS) {
              throw std::runtime_error("Failed to read input data");
            }
          }
          // Execute all algorithms
          for (auto& alg : m_algorithms) {
            StopWatch sw(localClocksAlgorithms[ialgo++]);
            if (alg->execute(++context) != ProcessCode::SUCCESS) {
              throw std::runtime_error("Failed to process event data");
            }
          }
          // Write out results
          for (auto& wrt : m_writers) {
            StopWatch sw(localClocksAlgorithms[ialgo++]);
            if (wrt->write(++context) != ProcessCode::SUCCESS) {
              throw std::runtime_error("Failed to write output data");
            }
          }
          ACTS_INFO("finished event " << event);
        }

        // add timing info to global information
        {
          tbb::queuing_mutex::scoped_lock lock(clocksAlgorithmsMutex);
          for (std::size_t i = 0; i < clocksAlgorithms.size(); ++i) {
            clocksAlgorithms[i] += localClocksAlgorithms[i];
          }
        }
      });

  // run end-of-run hooks
  for (auto& wrt : m_writers) {
    names.push_back("Writer:" + wrt->name() + ":endRun");
    clocksAlgorithms.push_back(Duration::zero());
    StopWatch sw(clocksAlgorithms.back());
    if (wrt->endRun() != ProcessCode::SUCCESS) { return EXIT_FAILURE; }
  }

  // summarize timing
  Duration totalWall = Clock::now() - clockWallStart;
  Duration totalReal = std::accumulate(
      clocksAlgorithms.begin(), clocksAlgorithms.end(), Duration::zero());
  std::size_t numEvents = endEvent - m_cfg.skip;
  ACTS_INFO("Processed " << numEvents << " events in " << asString(totalWall)
                         << " (wall clock)");
  ACTS_INFO("Average time per event: " << perEvent(totalReal, numEvents));
  ACTS_DEBUG("Average time per algorithm:");
  for (size_t i = 0; i < names.size(); ++i) {
    ACTS_DEBUG("  " << names[i] << ": "
                    << perEvent(clocksAlgorithms[i], numEvents));
  }
  storeTiming(names,
              clocksAlgorithms,
              numEvents,
              joinPaths(m_cfg.outputDir, "timing.tsv"));

  return EXIT_SUCCESS;
}
