#!/bin/bash
docker start test_fvib_casm_container
docker exec -it -w $PWD test_fvib_casm_container bash
