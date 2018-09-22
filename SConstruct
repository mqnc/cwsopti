
env = Environment(TARGET_ARCH = 'x86')
env.VariantDir("build", "src")
env.AppendUnique(CPPPATH=["src"])

if "msvc" in env["TOOLS"]:
	env.AppendUnique(CXXFLAGS=["/O2"])
	#env.AppendUnique(CXXFLAGS=["/DEBUG"])
	env.AppendUnique(CXXFLAGS=["/EHsc"])

else:
	#env.AppendUnique(CCFLAGS=["-Wall", "-Wextra", "-g", "-pedantic", "-pthread"])
	env.AppendUnique(CFLAGS=["-std=gnu99"])
	env.AppendUnique(CXXFLAGS=["-std=c++14"])
	env.AppendUnique(LINKFLAGS=["-pthread"])

env.Program("program", ["build/main.cpp", "build/sh/spherical_harmonics.cc", "build/sh/default_image.cc"])
