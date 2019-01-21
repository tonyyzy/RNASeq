# RNASeq

[![Build Status](https://travis-ci.com/tonyyzy/RNASeq.svg?token=5Fwptxoz1iaezXoMzRSd&branch=master)](https://travis-ci.com/tonyyzy/RNASeq)

## Testing
We use `python` and `pytest` as the testing framework. All the tests are stored in `./tests`. Before committing a new test, you should test your test locally first. 

Assuming you have all dependencies met

run `$ pytest ./tests`

You should see something similar to the following:

```
$ pytest tests
============================= test session starts ==============================
platform linux -- Python 3.6.3, pytest-3.3.0, py-1.5.2, pluggy-0.6.0
rootdir: /home/travis/build/tonyyzy/RNASeq, inifile:
collected 1 item                                                               
tests/test_STAR_index_nodocker.py .                                      [100%]
=========================== 1 passed in 2.01 seconds ===========================
The command "pytest tests" exited with 0.
```