#ifndef SPINNER_H
#define SPINNER_H

#include <iostream>
#include <thread>
#include <chrono>
#include <atomic>
#include <string>


/**
 * @brief Class for displaying a spinner in the console while a process is running
 * @note The class is not thread-safe, so only one instance of the class should be used at a time
 * @note It has the following functions:
 *      - Spinner(const std::string& processMsg, const std::string& timeMsg = "Total time: ", int intervalMs = 100)
 *          - processMsg: the message to display before the spinner
 *          - timeMsg: the message to display after the spinner
 *          - intervalMs: the interval between each spinner frame (in milliseconds)
 *      - ~Spinner()
 *          - Stop the spinner and print the total time
 */
class Spinner
{
public:
    using Clock = std::chrono::steady_clock;

    explicit Spinner(const std::string& processMsg, 
                     const std::string& timeMsg = "Total time: ", 
                     int intervalMs = 100)
        : running_(true),
          process_msg_(processMsg),
          time_msg_(timeMsg),
          interval_ms_(intervalMs),
          start_time_(Clock::now())
    {
       
        spinner_thread_ = std::thread([this]() {
            static const char spinner_chars[] = { '|', '/', '-', '\\' };
            int i = 0;
            while (running_)
            {
                
                std::cout << "\r" 
                          << process_msg_ << spinner_chars[i++] 
                          << std::flush;
                i %= 4;
                std::this_thread::sleep_for(std::chrono::milliseconds(interval_ms_));
            }
        });
    }
    ~Spinner()
    {
        stop();
    }

    void stop()
    {
        if (running_)
        {
            running_ = false;
            if (spinner_thread_.joinable())
            {
                spinner_thread_.join();
            }
            std::cout << "\r\033[2K\r" << std::flush;  
            auto end_time = Clock::now();
            double elapsed = std::chrono::duration<double>(end_time - start_time_).count();
            std::cout << "\033[1;32m[Time]\033[0m " 
                      << time_msg_ 
                      << "\033[1;32m" << elapsed << "s\033[0m" 
                      << std::endl;
        }
    }

private:
    std::atomic<bool> running_;
    std::thread spinner_thread_;

    std::string process_msg_;  // Display message before the spinner, e.g. "Processing..."
    std::string time_msg_;     // Time message to display after the spinner, e.g. "Total time: "
    int interval_ms_;          // refresh interval

    std::chrono::time_point<Clock> start_time_; // Record the start time
};

#endif // SPINNER_H
