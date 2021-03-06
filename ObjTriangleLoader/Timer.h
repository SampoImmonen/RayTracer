#pragma once
#include <chrono>
#include <iostream>

/*
A simple scope based timer. Currently only a header file
To be expanded
*/


class Timer {

public:
	Timer(bool toFile = false)
	{
		writetoFile = toFile;
		m_StartTimepoint = std::chrono::high_resolution_clock::now();
	}

	~Timer()
	{
		Stop();
	}

private:



	void Stop() {
		auto EndTimepoint = std::chrono::high_resolution_clock::now();

		auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimepoint).time_since_epoch().count();
		auto end = std::chrono::time_point_cast<std::chrono::microseconds>(EndTimepoint).time_since_epoch().count();

		auto duration = end - start;
		double ms = duration * 0.001;
		std::cout << "Timer:\n";
		std::cout << " " << duration << "us\n";
		std::cout << " " << ms << "ms\n";
	}
	std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTimepoint;
	bool writetoFile;
};

