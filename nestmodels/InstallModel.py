"""Install adaptive LIF neuron using nestml"""

import nest
from pynestml.frontend.pynestml_frontend import to_nest, install_nest
NEST_SIMULATOR_INSTALL_LOCATION = nest.ll_api.sli_func("statusdict/prefix ::")
to_nest(input_path="adapt_lif.nestml",
        target_path="/tmp/nestml-component",
        logging_level="ERROR")
install_nest("/tmp/nestml-component", NEST_SIMULATOR_INSTALL_LOCATION)
