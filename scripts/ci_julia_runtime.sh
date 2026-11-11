#!/usr/bin/env bash
set -euo pipefail

artifact=$2
case "$1" in
  pack)
    root=$(julia --startup-file=no -e '
      @assert read(joinpath(ARGS[1], "julia-version.txt"), String) == string(VERSION, " ", Base.GIT_VERSION_INFO.commit)
      print(dirname(Sys.BINDIR))' "$artifact")
    [[ "$root" == "$RUNNER_TOOL_CACHE/"* ]]
    printf '%s\n' "${root#"$RUNNER_TOOL_CACHE/"}" > "$artifact/julia-runtime-path.txt"
    tar --exclude='./share/doc' --exclude='./share/man' -cf - -C "$root" . |
      gzip -1 > "$artifact/julia-runtime.tar.gz"
    ;;
  restore)
    IFS= read -r relative < "$artifact/julia-runtime-path.txt"
    case "$relative" in ''|/*|*..*) exit 1 ;; esac
    root="$RUNNER_TOOL_CACHE/$relative"
    mkdir -p "$root"
    tar -xzf "$artifact/julia-runtime.tar.gz" -C "$root"
    printf '%s\n' "$root/bin" >> "$GITHUB_PATH"
    ;;
  *) exit 1 ;;
esac
