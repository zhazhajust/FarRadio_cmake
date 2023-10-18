## Prerequisites

* A compiler with C++11 support
* Pip 10+ or CMake >= 3.4 (or 3.14+ on Windows, which was the first version to support VS 2019)
* Ninja or Pip 10+


## Installation

Just clone this repository and pip install. Note the `--recursive` option which is
needed for the pybind11 submodule:

```bash
git clone --recursive https://github.com/zhazhajust/Faradio_cmake.git
python setup.py install
```

With the `setup.py` file included in this example, the `pip install` command will
invoke CMake and build the pybind11 module as specified in `CMakeLists.txt`.


## License

Faradio_cmake is provided under a BSD-style license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.


## Test call

```python
import pyfaradio
assert pyfaradio.__version__ == "0.0.1"
```

./pyfaradio.faradio.run
