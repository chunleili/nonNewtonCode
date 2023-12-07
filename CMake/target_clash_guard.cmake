macro(target_clash_guard target_name target_dir)
    if(NOT TARGET ${target_name})
        add_subdirectory(${target_dir})
    else()
        message(STATUS "${target_name} already exists, use the existing one")
    endif()
endmacro(target_clash_guard)