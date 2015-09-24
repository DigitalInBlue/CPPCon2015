#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <bitset>
#include <cmath>

#ifndef WIN32
	#define WIN32 0
#endif

/// Use to find out the binary representation of a float.
typedef union 
{
	float number;
  
	struct 
	{
		uint32_t mantisa : 23;
		uint32_t exponent : 8;
		uint32_t sign : 1;
	} parts;

	uint32_t raw;
} ieeeFloat;

/// Use to find out the binary representation of a double.
typedef union 
{
	double number;
  
	struct 
	{
		uint64_t mantisa : 52;
		uint64_t exponent : 11;
		uint64_t sign : 1;
	} parts;

	uint64_t raw;
} ieeeDouble;

/// Format nice output for google test.
std::ostream& cppCon(const int widthMult = 2)
{
	std::cerr << std::left << "[  CPPCON  ] " << ::testing::UnitTest::GetInstance()->current_test_info()->test_case_name() << ".";
	std::cerr << ::testing::UnitTest::GetInstance()->current_test_info()->name() << " -- ";
	std::cerr << std::fixed << std::right;
	std::cerr << std::setw(std::numeric_limits<double>::digits10 * widthMult) << std::setprecision(std::numeric_limits<double>::digits10 * widthMult);
	return std::cerr;
}

/// Heron's formula to compute the area of a triangle.
template<typename T> T HeronsFormula(const T a, const T b, const T c)
{
	const auto s = ((a + b) + c) / 2;
	cppCon() << "S = " << s << "\n";
	const auto num = s * (s - a) * (s - b) * (s - c);
	const auto area = std::sqrt(num);
	return area;
}

/// Kahan's formulat to compute the area of a triangle
template<typename T> T KahansFormula(const T a, const T b, const T c)
{
	assert((a >= b) && (b >= c));
	const auto num = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));
	const auto area = std::sqrt(num) / T(4);
	return area;
}

/// How is Pi stored?
TEST(CPPCon2015, pi)
{
	ieeeFloat floatPi;
	floatPi.number = 3.14159265358979323846264338327950288419716939937510582f;

	ieeeDouble doublePi;
	doublePi.number = 3.14159265358979323846264338327950288419716939937510582;

	cppCon() << "Float Pi: " << floatPi.number << "\n";
	{
		std::bitset<1> sign = floatPi.parts.sign;
		std::bitset<8> exponent = floatPi.parts.exponent;
		std::bitset<23> mantisa = floatPi.parts.mantisa;
		cppCon() << "Float Pi: [" << sign << "][" << exponent << "][" << mantisa << "]\n";
	}

	cppCon() << "Double Pi: " << doublePi.number << "\n";
	{
		std::bitset<1> sign = doublePi.parts.sign;
		std::bitset<11> exponent = doublePi.parts.exponent;
		std::bitset<52> mantisa = doublePi.parts.mantisa;
		cppCon() << "Double Pi: [" << sign << "][" << exponent << "][" << mantisa << "]\n";
	}
}

/// How is Pi stored?
#if(WIN32 && _MSC_VER > 1700) 
TEST(CPPCon2015, piNext)
{
	ieeeFloat floatPi;
	floatPi.number = 3.14159265358979323846264338327950288419716939937510582f;

	cppCon() << "Float Pi: " << floatPi.number << "\n";

	{
		std::bitset<1> sign = floatPi.parts.sign;
		std::bitset<8> exponent = floatPi.parts.exponent;
		std::bitset<23> mantisa = floatPi.parts.mantisa;
		cppCon() << "Float Pi: [" << sign << "][" << exponent << "][" << mantisa << "]\n";
	}
	
	// Get the next float after PI in the direction of 3.0.
	floatPi.number = std::nextafterf(floatPi.number, 3.0);

	cppCon() << "Float Pi Next: " << floatPi.number << "\n";

	{
		std::bitset<1> sign = floatPi.parts.sign;
		std::bitset<8> exponent = floatPi.parts.exponent;
		std::bitset<23> mantisa = floatPi.parts.mantisa;
		cppCon() << "Float Pi Next: [" << sign << "][" << exponent << "][" << mantisa << "]\n";
	}
}
#endif

/// 0.1 + 0.2 == 0.3...never
TEST(CPPCon2015, pointOnePlusPointTwo)
{
	auto zeroPointOne = 0.1f;
	auto zeroPointTwo = 0.2f;
	auto zeroPointThree = 0.3f;
	auto sum = zeroPointOne + zeroPointTwo;

	EXPECT_DOUBLE_EQ(0.3f, zeroPointThree);
	EXPECT_EQ(0.3f, zeroPointThree);
	EXPECT_EQ(0.3f, sum);
	cppCon() << "zeroPointOne == " << zeroPointOne << "\n";
	cppCon() << "zeroPointTwo == " << zeroPointTwo << "\n";
	cppCon() << "zeroPointThree == " << zeroPointThree << "\n";
	cppCon() << "sum == " << sum << "\n";
	EXPECT_EQ(zeroPointThree, sum);
	EXPECT_DOUBLE_EQ(zeroPointThree, sum);
}

/// 0.1 + 0.2 == 0.3...never
TEST(CPPCon2015, pointOnePlusPointTwoDoubles)
{
	auto zeroPointOne = 0.1;
	auto zeroPointTwo = 0.2;
	auto zeroPointThree = 0.3;
	auto sum = zeroPointOne + zeroPointTwo;

	EXPECT_DOUBLE_EQ(0.3, zeroPointThree);
	EXPECT_EQ(0.3, zeroPointThree);
	EXPECT_EQ(0.3, sum) << "This may pass or may fail depending on your floating point switch when compiling.";

	cppCon() << "zeroPointOne == " << zeroPointOne << "\n";
	cppCon() << "zeroPointTwo == " << zeroPointTwo << "\n";
	cppCon() << "zeroPointThree == " << zeroPointThree << "\n";
	cppCon() << "sum == " << sum << "\n";
	EXPECT_NE(zeroPointThree, sum);
	EXPECT_DOUBLE_EQ(zeroPointThree, sum);
}

/// How is 1.0f stored?
TEST(CPPCon2015, one)
{
	ieeeFloat float1;
	float1.number = 1.0f;
	
	cppCon() << "Float 1.0:  " << float1.number << "\n";

	{
		std::bitset<1> sign = float1.parts.sign;
		std::bitset<8> exponent = float1.parts.exponent;
		std::bitset<23> mantisa = float1.parts.mantisa;
		cppCon() << "Float 1.0: [" << sign << "][" << exponent << "][" << mantisa << "]\n";
	}

	EXPECT_DOUBLE_EQ(1.0, float1.number);

	float1.parts.mantisa += 1;

	{
		std::bitset<1> sign = float1.parts.sign;
		std::bitset<8> exponent = float1.parts.exponent;
		std::bitset<23> mantisa = float1.parts.mantisa;
		cppCon() << "Float 1.0 + \"1\": [" << sign << "][" << exponent << "][" << mantisa << "]\n";
	}

	cppCon() << "Float 1.0 + \"1\":  " << float1.number << "\n";
}

/// How is the very small stored?
TEST(CPPCon2015, verySmalll)
{
	EXPECT_EQ(-37, FLT_MIN_10_EXP);

	ieeeFloat float1;
	float1.number = 1.0e-38f;

	cppCon(3) << "Float 1.0e-38:  " << float1.number << "\n";

	{
		std::bitset<1> sign = float1.parts.sign;
		std::bitset<8> exponent = float1.parts.exponent;
		std::bitset<23> mantisa = float1.parts.mantisa;
		cppCon() << "Float 1.0e-38: [" << sign << "][" << exponent << "][" << mantisa << "]\n";
	}
}

#if(WIN32 && _MSC_VER > 1700) 
/// How is 0.1f stored?
TEST(CPPCon2015, special)
{
	ieeeFloat float1;

	// Divide by zero
	{
		float x = 1.0f;
		float y = 1.0f;
		float1.number = (1.0f) / (x - y);
		std::bitset<1> sign = float1.parts.sign;
		std::bitset<8> exponent = float1.parts.exponent;
		std::bitset<23> mantisa = float1.parts.mantisa;
		cppCon() << "DIV 0: [" << sign << "][" << exponent << "][" << mantisa << "]\n";
	}

	// NaN
	{
		float1.number = NAN;
		std::bitset<1> sign = float1.parts.sign;
		std::bitset<8> exponent = float1.parts.exponent;
		std::bitset<23> mantisa = float1.parts.mantisa;
		cppCon() << "NAN: [" << sign << "][" << exponent << "][" << mantisa << "]\n";
	}

	// Infinity
	{
		float1.number = INFINITY;
		std::bitset<1> sign = float1.parts.sign;
		std::bitset<8> exponent = float1.parts.exponent;
		std::bitset<23> mantisa = float1.parts.mantisa;
		cppCon() << "INFINITY: [" << sign << "][" << exponent << "][" << mantisa << "]\n";
	}
}
#endif

/// Simulation time, accumulating at 100 Hz
TEST(CPPCon2015, simTime100HzFloatAccumulate)
{
	auto totalFrames = size_t(0);
	auto frameLength = 0.01f;
	auto simTime = 0.0f;

	// Run 120 Frames
	auto simTime120Frames = 1.20f;
	for(; totalFrames < 120; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_NE(simTime120Frames, simTime);
	cppCon() << simTime120Frames << " != " << simTime << " (" << simTime120Frames - simTime << ")\n";

	// Run 1200 Frames
	auto simTime1200Frames = 12.00f;
	for(; totalFrames < 1200; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_NE(simTime1200Frames, simTime);
	cppCon() << simTime1200Frames << " != " << simTime << " (" << simTime1200Frames - simTime << ")\n";

	// Run 12000 Frames
	auto simTime12000Frames = 120.0f;
	for(; totalFrames < 12000; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_NE(simTime12000Frames, simTime);
	cppCon() << simTime12000Frames << " != " << simTime << " (" << simTime12000Frames - simTime << ")\n";

	// Run 60000 Frames
	auto simTime60000Frames = 600.0f;
	for(; totalFrames < 60000; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_NE(simTime60000Frames, simTime);
	cppCon() << simTime60000Frames << " != " << simTime << " (" << simTime60000Frames - simTime << ")\n";

	// Run 360000 Frames
	auto simTime360000Frames = 3600.0f;
	for(; totalFrames < 360000; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_NE(simTime360000Frames, simTime);
	cppCon() << simTime360000Frames << " != " << simTime << " (" << simTime360000Frames - simTime << ")\n";
}

/// Simulation time, accumulating at 100 Hz
TEST(CPPCon2015, simTime100HzFloatDelta)
{
	auto totalFrames = size_t(0);
	auto frameLength = 0.01f;
	auto simTime = 0.0f;

	// Run 120 Frames
	auto simTime120Frames = 1.20f;
	totalFrames = 120;
	simTime = totalFrames * frameLength;
	EXPECT_NE(simTime120Frames, simTime);
	cppCon() << simTime120Frames << " != " << simTime << " (" << simTime120Frames - simTime << ")\n";

	// Run 1200 Frames
	auto simTime1200Frames = 12.00f;
	totalFrames = 1200;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime1200Frames, simTime);
	cppCon() << simTime1200Frames << " == " << simTime << " (" << simTime1200Frames - simTime << ")\n";

	// Run 12000 Frames
	auto simTime12000Frames = 120.0f;
	totalFrames = 12000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime12000Frames, simTime);
	cppCon() << simTime12000Frames << " == " << simTime << " (" << simTime12000Frames - simTime << ")\n";

	// Run 60000 Frames
	auto simTime60000Frames = 600.0f;
	totalFrames = 60000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime60000Frames, simTime);
	cppCon() << simTime60000Frames << " == " << simTime << " (" << simTime60000Frames - simTime << ")\n";

	// Run 360000 Frames
	auto simTime360000Frames = 3600.0f;
	totalFrames = 360000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime360000Frames, simTime);
	cppCon() << simTime360000Frames << " == " << simTime << " (" << simTime360000Frames - simTime << ")\n";
}

/// Simulation time, accumulating at 100 Hz
TEST(CPPCon2015, simTime100HzDoubleAccumulate)
{
	auto totalFrames = size_t(0);
	auto frameLength = 0.01;
	auto simTime = 0.0;

	// Run 120 Frames
	auto simTime120Frames = 1.20;
	for(; totalFrames < 120; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_NE(simTime120Frames, simTime);
	cppCon() << simTime120Frames << " != " << simTime << " (" << simTime120Frames - simTime << ")\n";

	// Run 1200 Frames
	auto simTime1200Frames = 12.00;
	for(; totalFrames < 1200; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_NE(simTime1200Frames, simTime);
	cppCon() << simTime1200Frames << " != " << simTime << " (" << simTime1200Frames - simTime << ")\n";

	// Run 12000 Frames
	auto simTime12000Frames = 120.0;
	for(; totalFrames < 12000; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_NE(simTime12000Frames, simTime);
	cppCon() << simTime12000Frames << " != " << simTime << " (" << simTime12000Frames - simTime << ")\n";

	// Run 60000 Frames
	auto simTime60000Frames = 600.0;
	for(; totalFrames < 60000; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_NE(simTime60000Frames, simTime);
	cppCon() << simTime60000Frames << " != " << simTime << " (" << simTime60000Frames - simTime << ")\n";

	// Run 360000 Frames
	auto simTime360000Frames = 3600.0;
	for(; totalFrames < 360000; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_NE(simTime360000Frames, simTime);
	cppCon() << simTime360000Frames << " != " << simTime << " (" << simTime360000Frames - simTime << ")\n";
}

/// Simulation time, accumulating at 100 Hz
TEST(CPPCon2015, simTime100HzDoubleDelta)
{
	auto totalFrames = size_t(0);
	auto frameLength = 0.01;
	auto simTime = 0.0;

	// Run 120 Frames
	auto simTime120Frames = 1.20;
	totalFrames = 120;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime120Frames, simTime);

	// Run 1200 Frames
	auto simTime1200Frames = 12.00;
	totalFrames = 1200;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime1200Frames, simTime);

	// Run 12000 Frames
	auto simTime12000Frames = 120.0;
	totalFrames = 12000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime12000Frames, simTime);

	// Run 600000 Frames
	auto simTime60000Frames = 600.0;
	totalFrames = 60000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime60000Frames, simTime);
}

/// Simulation time, accumulating at 100 Hz
TEST(CPPCon2015, simTime100HzDoubleDeltaBig)
{
	auto totalFrames = size_t(0);
	auto frameLength = 0.01;
	auto simTime = 0.0;

	// Run 120 Frames
	auto simTime120Frames = 1.20;
	totalFrames = 120;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime120Frames, simTime);

	// Run 1200 Frames
	auto simTime1200Frames = 12.00;
	totalFrames = 1200;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime1200Frames, simTime);

	// Run 12000 Frames
	auto simTime12000Frames = 120.0;
	totalFrames = 12000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime12000Frames, simTime);

	auto simTime1200000Frames = 12000.0;
	totalFrames = 1200000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime1200000Frames, simTime);

	auto simTime120000000Frames = 1200000.0;
	totalFrames = 120000000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime120000000Frames, simTime);

	auto simTime120000000000Frames = 1200000000.0;
	totalFrames = 120000000000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime120000000000Frames, simTime);

	auto simTime120000000000000Frames = 1200000000000.0;
	totalFrames = 120000000000000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime120000000000000Frames, simTime);

	auto simTime120000000000000000Frames = 1200000000000000.0;
	totalFrames = 120000000000000000;
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime120000000000000000Frames, simTime);
}

/// A real example from a real simulation
TEST(CPPCon2015, simTimeFromTheWild)
{
	double simTimeDt = 0;
	double simTime = 0;

	// We have some model that runs at 100Hz.
	double modelLastUpdateTime = 0;
	volatile double modelUpdateDt = 100.0;

	size_t frameCount = 0;

	auto frameLoop = [&simTimeDt, &simTime, &modelLastUpdateTime, &modelUpdateDt, &frameCount](size_t frameCountLimit) 
		{
			while(frameCount < frameCountLimit)
			{
				// "query" our model
				simTimeDt = (modelLastUpdateTime + (1.0/modelUpdateDt)) - simTime;

				// Accumulate sim time.
				simTime += simTimeDt;

				// "update" our model
				modelLastUpdateTime = simTime;

				frameCount++;
			}
	};

	// Run 120 Frames
	auto simTime120Frames = 1.20;
	frameLoop(120);
	EXPECT_NE(simTime120Frames, simTime);
	cppCon() << simTime120Frames << " != " << simTime << " (" << simTime120Frames - simTime << ")\n";

	// Run 1200 Frames
	auto simTime1200Frames = 12.00;
	frameLoop(1200);
	EXPECT_NE(simTime1200Frames, simTime);
	cppCon() << simTime1200Frames << " != " << simTime << " (" << simTime1200Frames - simTime << ")\n";

	// Run 12000 Frames
	auto simTime12000Frames = 120.0;
	frameLoop(12000);
	EXPECT_NE(simTime12000Frames, simTime);
	cppCon() << simTime12000Frames << " != " << simTime << " (" << simTime12000Frames - simTime << ")\n";

	// Run 600000 Frames
	auto simTime60000Frames = 600.0;
	frameLoop(60000);
	EXPECT_NE(simTime60000Frames, simTime);
	cppCon() << simTime60000Frames << " != " << simTime << " (" << simTime60000Frames - simTime << ")\n";
}

/// Simulation time, accumulating at sub-microseconds using floats.
TEST(CPPCon2015, simTimeSubMicrosecondFloatDelta)
{
	auto totalFrames = 0e0f;
	auto frameLength = 1.2345e-6f;
	auto simTime = 0.0f;

	for(; totalFrames < 12000; totalFrames++)
	{
		simTime += frameLength;
	}

	// Run 10.0e3 Frames
	auto simTime10e3Frames = 1.2345e-2f;
	for(; totalFrames < 10e3f; totalFrames++)
	{
		simTime += frameLength;
	}
	simTime = totalFrames * frameLength;
	EXPECT_EQ(simTime10e3Frames, simTime);
	cppCon() << simTime10e3Frames << " != " << simTime << " (" << simTime10e3Frames - simTime << ")\n";

	// Run 10.0e6 Frames
	auto simTime10e6Frames = 12.345f;
	for(; totalFrames < 10e6; totalFrames++)
	{
		simTime += frameLength;
	}
	EXPECT_EQ(simTime10e6Frames, simTime);
	cppCon() << simTime10e6Frames << " != " << simTime << " (" << simTime10e6Frames - simTime << ")\n";
}

TEST(CPPCon2015, Inverse1)
{
	auto x = 299792458.5f;
	EXPECT_NE(x / 10.0f, x * 0.1f);

	EXPECT_NE(29979245.85f, x / 10.0f);

	// Google Test didn't even catch this.  It thought things were fine. // EXPECT_NE(29979245.85f, x * 0.1f);
	EXPECT_NE(29979245.85, static_cast<double>(x * 0.1f));

	cppCon() << "x / 10.0f = " << x / 10.0f << " (" << 29979245.85f - (x / 10.0f) << ")\n";
	cppCon() << "x * 0.1f  = " << x * 0.1f << " (" << 29979245.85 - static_cast<double>(x * 0.1f) << ")\n";
}

TEST(CPPCon2015, Inverse2)
{
	auto x = 299792458.5f;
	EXPECT_NE(x / 10.0f, x * 0.1f);

	EXPECT_NE(2997924585.0f, x * 10.0f);
	EXPECT_NE(2997924585.0f, x / 0.1f);

	cppCon() << "x * 10.0f = " << x * 10.0f << " (" << 2997924585.0f - (x * 10.0f) << ")\n";
	cppCon() << "x / 0.1f  = " << x / 0.1f << " (" << 2997924585.0f - (x / 0.1f) << ")\n";
}

/// Area of a triangle
/// From http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
/// From https://www.cs.berkeley.edu/~wkahan/Triangle.pdf
TEST(CPPCon2015, AreaOfATriangleFloat)
{
	const auto a = 100000.0f;
	const auto b = 99999.99979f;
	const auto c = 0.00029f;

	ASSERT_TRUE(a >= b);
	ASSERT_TRUE(b >= c);

	auto heronsFormula = HeronsFormula(a, b, c);
	auto kahansFormula = KahansFormula(a, b, c);
	EXPECT_NE(kahansFormula, heronsFormula);
	cppCon() << "Heron: " << heronsFormula << "\n";
	cppCon() << "Kahan: " << kahansFormula << "\n";
	cppCon() << "delta: " << kahansFormula - heronsFormula << "\n";
}

/// Area of a triangle
/// From http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
/// From https://www.cs.berkeley.edu/~wkahan/Triangle.pdf
TEST(CPPCon2015, AreaOfATriangleDouble)
{
	const auto a = 100000.0;
	const auto b = 99999.99979;
	const auto c = 0.00029;

	ASSERT_TRUE(a >= b);
	ASSERT_TRUE(b >= c);

	auto heronsFormula = HeronsFormula(a, b, c);
	auto kahansFormula = KahansFormula(a, b, c);
	EXPECT_NE(kahansFormula, heronsFormula);
	cppCon() << "Heron: " << heronsFormula << "\n";
	cppCon() << "Kahan: " << kahansFormula << "\n";
	cppCon() << "delta: " << kahansFormula - heronsFormula << "\n";
}

/// xkcd thought it was funny.
TEST(CPPCon2015, xkcd)
{
	std::setprecision(32);
	
	const auto pi = 3.14159265359;
	const auto e = 2.71828182845904;

	EXPECT_NE(20.0, std::pow(e, pi) - pi);
	EXPECT_NEAR(20.0, std::pow(e, pi) - pi, 0.001);
}

#if(WIN32 && _MSC_VER > 1700) 
TEST(CppCon2015, next)
{
	auto max = 10.0f;
	size_t counter = 0;
	auto start = 0.0f;

	for(float f = start; f < 10.0f; f = std::nextafter(f, max))
	{
		counter++;

		if(std::fmodf(f, 0.1f) == 0)
		{
			cppCon() << "Number of floats from " << start << " to " << f << " = " << counter << "\n";
			counter = 0;
			start = f;
		}
	}
}
#endif

TEST(CppCon2015, radian)
{
	const auto oneRadian = 0.15915494309f;
	const auto control = 0.000000000015915494309f;

	const auto oneRadianMultiplied = oneRadian * 1.0e-10f;
	const auto oneRadianDivided = oneRadian / 1.0e10f;

	cppCon() << "Control: " << control << "\n";
	cppCon() << "x*1.0e-10f: " << oneRadianMultiplied << "(" << fabs(oneRadianMultiplied - control) << ")\n";
	cppCon() << "x/1.0e10f: " << oneRadianDivided << "(" << fabs(oneRadianDivided - control) << ")\n";
	cppCon() << "Relative Error: " << fabs(oneRadianDivided - control) / control << "\n";

	EXPECT_EQ(oneRadianDivided, oneRadianMultiplied);
}
