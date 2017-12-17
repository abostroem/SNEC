#!/usr/bin/env python3

import check_run
import prep_snec

parameters = prep_snec.get_input_parameters()
search_dir = '/home/bostroem/SNEC/snec_models'

check_run.check_run(parameters, search_dir)
