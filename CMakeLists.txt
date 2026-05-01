if(CONFIG_CTL_TIMER)
	target_include_directories(app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/controller/timer)
endif()

if(CONFIG_FLT_LPF)
	target_include_directories(app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/filter/lpf)
endif()

if(CONFIG_ID_RLS)
	target_sources(app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/identify/rls/rls.cpp)
	target_include_directories(app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/identify/rls)
endif()

if(CONFIG_CTL_PID)
	target_sources(app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/controller/pid/pid.cpp)
	target_include_directories(app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/controller/pid)
endif()

if(CONFIG_CTL_POWER)
	target_sources(app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/controller/power_ctrl/power_ctrl.cpp)
	target_include_directories(app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/controller/power_ctrl)
endif()



