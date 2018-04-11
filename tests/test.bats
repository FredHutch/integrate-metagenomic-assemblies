#!/usr/bin/env bats

@test "aws-batch-helpers code" {
  python3 -c 'from batch_helpers.helpers import run_cmds; run_cmds(["echo", "success"], stdout="/usr/local/tests/installed")'
  
  [[ "$(cat /usr/local/tests/installed)" == "success" ]]
}
