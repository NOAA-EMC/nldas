baseforcing_pluginMod.o baseforcing_pluginMod.d : baseforcing_pluginMod.F90
baseforcing_pluginMod.o : geosdomain_module.o
baseforcing_pluginMod.o : gdasdomain_module.o
baseforcing_pluginMod.o : ecmwfdomain_module.o
baseforcing_pluginMod.o : nldasdomain_module.o
baseforcing_pluginMod.o : nldas2domain_module.o
baseforcing_pluginMod.o : gswpdomain_module.o
baseforcing_pluginMod.o : bergdomain_module.o
