/* Copyright 2022 Philip Miller
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "SignalHandling.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_IParser.H>

#include <cctype>

// For sigaction() et al.
#if defined(__linux__) || defined(__APPLE__)
#include <signal.h>
#endif

std::atomic<bool> SignalState::signal_received_flags[NUM_SIGNALS];
bool SignalState::signal_conf_requests_break[NUM_SIGNALS];
bool SignalState::signal_conf_requests_checkpoint[NUM_SIGNALS];
bool SignalState::signal_actions_requested[2];
#if defined(AMREX_USE_MPI)
MPI_Request SignalState::signal_mpi_ibcast_request;
#endif

int
SignalState::parseSignalNameToNumber(const std::string &str)
{
    amrex::IParser signals_parser(str);

#if defined(__linux__) || defined(__APPLE__)
    struct {
        const char* abbrev;
        const int value;
    } signals_to_parse[] = {
        {"ABRT", SIGABRT},
        {"ALRM", SIGALRM},
        {"BUS", SIGBUS},
        {"CHLD", SIGCHLD},
        {"CLD", SIGCHLD}, // Synonymous to SIGCHLD on Linux
        {"CONT", SIGCONT},
#if defined(SIGEMT)
        {"EMT", SIGEMT}, // macOS and some Linux architectures
#endif
        // Omitted because AMReX typically handles SIGFPE specially
        // {"FPE", SIGFPE},
        {"HUP", SIGHUP},
        {"ILL", SIGILL},
#if defined(SIGINFO)
        {"INFO", SIGINFO}, // macOS and some Linux architectures
#endif
        {"INT", SIGINT},
        {"IO", SIGIO},
        {"IOT", SIGABRT}, // Synonymous to SIGABRT on Linux
        // {"KILL", SIGKILL}, // Cannot be handled
        {"PIPE", SIGPIPE},
        {"POLL", SIGIO}, // Synonymous to SIGIO on Linux
        {"PROF", SIGPROF},
#if defined(SIGPWR)
        {"PWR", SIGPWR}, // Linux-only
#endif
        {"QUIT", SIGQUIT},
        {"SEGV", SIGSEGV},
#if defined(SIGSTKFLT)
        {"STKFLT", SIGSTKFLT}, // Linux-only
#endif
        // {"STOP", SIGSTOP}, // Cannot be handled
        {"SYS", SIGSYS},
        {"TERM", SIGTERM},
        {"TRAP", SIGTRAP},
        {"TSTP", SIGTSTP},
        {"TTIN", SIGTTIN},
        {"TTOU", SIGTTOU},
        {"URG", SIGURG},
        {"USR1", SIGUSR1},
        {"USR2", SIGUSR2},
        {"VTALRM", SIGVTALRM},
        {"WINCH", SIGWINCH},
        {"XCPU", SIGXCPU},
        {"XFSZ", SIGXFSZ},
    };

    for (const auto& sp : signals_to_parse) {
        std::string name_upper = sp.abbrev;
        std::string name_lower = name_upper;
        for (char &c : name_lower) {
            c = std::tolower(c);
        }

        signals_parser.setConstant(name_upper, sp.value);
        signals_parser.setConstant(name_lower, sp.value);
        name_upper = "SIG" + name_upper;
        name_lower = "sig" + name_lower;
        signals_parser.setConstant(name_upper, sp.value);
        signals_parser.setConstant(name_lower, sp.value);
    }
#endif // #if defined(__linux__) || defined(__APPLE__)

    auto spf = signals_parser.compileHost<0>();

    int sig = spf();
    AMREX_ALWAYS_ASSERT(sig < NUM_SIGNALS);

    return sig;
}

void
SignalState::InitSignalHandling()
{
#if defined(__linux__) || defined(__APPLE__)
    struct sigaction sa;
    sigemptyset(&sa.sa_mask);
    for (int i = 0; i < NUM_SIGNALS; ++i) {
        signal_received_flags[i] = false;
        if (signal_conf_requests_checkpoint[i] || signal_conf_requests_break[i]) {
            if (amrex::ParallelDescriptor::MyProc() == 0) {
                sa.sa_handler = &SignalState::SignalSetFlag;
            } else {
                sa.sa_handler = SIG_IGN;
            }
            int result = sigaction(i, &sa, nullptr);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(result == 0,
                                             "Failed to install signal handler for a configured signal");
        }
    }
#endif
}

void
SignalState::CheckSignals()
{
    // We assume that signals will definitely be delivered to rank 0,
    // and may be delivered to other ranks as well. For coordination,
    // we process them according to when they're received by rank 0.
    if (amrex::ParallelDescriptor::MyProc() == 0) {
        for (int i = 0; i < NUM_SIGNALS; ++i) {
            // Read into a local temporary to ensure the same value is
            // used throughout. Atomically exchange it with false to
            // unset the flag without risking loss of a signal - if a
            // signal arrives after this, it will be handled the next
            // time this function is called.
            bool signal_i_received = signal_received_flags[i].exchange(false);

            if (signal_i_received) {
                signal_actions_requested[SIGNAL_REQUESTS_BREAK] |= signal_conf_requests_break[i];
                signal_actions_requested[SIGNAL_REQUESTS_CHECKPOINT] |= signal_conf_requests_checkpoint[i];
            }
        }
    }

#if defined(AMREX_USE_MPI)
    auto comm = amrex::ParallelDescriptor::Communicator();
    MPI_Ibcast(signal_actions_requested, SIGNAL_REQUESTS_MAX+1, MPI_CXX_BOOL, 0, comm, &signal_mpi_ibcast_request);
#endif
}

void
SignalState::WaitSignals()
{
#if defined(AMREX_USE_MPI)
    MPI_Wait(&signal_mpi_ibcast_request, MPI_STATUS_IGNORE);
#endif
}

bool
SignalState::TestAndResetActionRequestFlag(int action_to_test)
{
    bool retval = signal_actions_requested[action_to_test];
    signal_actions_requested[action_to_test] = false;
    return retval;
}

void
SignalState::SignalSetFlag(int signal_number)
{
    signal_received_flags[signal_number] = true;
}

