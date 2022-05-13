add_executable(main main.cc)
target_link_libraries(main PUBLIC getInput cpimaging ROOT::Gpad)
set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")