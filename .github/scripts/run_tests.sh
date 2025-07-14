#!/bin/bash

output=$(./build/genesis4)
if [[ "$output" != *"Usage: genesis4"* ]]; then
  echo "Error: Expected 'Usage: genesis4' in output but got: $output"
  exit 1
else
  echo "Genesis4 seems to have built correctly; 'usage' found in output."
fi

# TODO: insert some actual tests here
