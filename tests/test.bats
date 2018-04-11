#!/usr/bin/env bats

@test "aws-batch-helpers code" {
  python3 -c 'from batch_helpers.helpers import run_cmds; run_cmds(["echo", "success"], stdout="/usr/local/ima/tests/installed")'
  
  [[ "$(cat /usr/local/ima/tests/installed)" == "success" ]]
}

@test "AWS CLI" {
  v="$(aws --version 2>&1)"
  [[ "$v" =~ "aws-cli" ]]
}

@test "Curl v7.47.0" {
  v="$(curl --version)"
  [[ "$v" =~ "7.47.0" ]]
}

@test "DIAMOND v0.9.10" {
  v="$(diamond --version)"
  [[ "$v" =~ "0.9.10" ]]
}

@test "MMseqs2 Release 2-23394" {
  v="$(mmseqs)"
  [[ "$v" =~ "2339462c06eab0bee64e4fc0ebebf7707f6e53fd" ]]
}

@test "Script in path" {
  h="$(integrate_assemblies.py -h || true)"
  echo "$h"
  [[ "$h" =~ "assemblies" ]]
}