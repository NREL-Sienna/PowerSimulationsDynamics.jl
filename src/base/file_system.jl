function check_folder(folder::String)
    !isdir(folder) && throw(IS.ConflictingInputsError("Specified folder is not valid"))
    try
        fake_dir =
            length(ARGS) > 0 ? joinpath(folder, string("fake_", ARGS[1])) :
            joinpath(folder, "fake")
        mkdir(fake_dir)
        rm(fake_dir)
    catch e
        throw(IS.ConflictingInputsError("Specified folder does not have write access [$e]"))
    end
end
