# Function to silently update the PID in launch.json
update_launch_json_pid = function() {
  current_dir = getwd()
  launch_json_path = file.path(current_dir, ".vscode", "launch.json")

  if (!file.exists(launch_json_path)) {
    return()
  }
  
  current_pid = Sys.getpid()
  json_content = readLines(launch_json_path, warn = FALSE)
  json_content = gsub('"pid": .*', paste0('"pid": ', current_pid, ","), json_content)
  writeLines(json_content, launch_json_path)
}

update_launch_json_pid()