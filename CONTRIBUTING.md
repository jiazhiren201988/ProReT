# Contributing to ProReT

Bug reports and focused feature requests are welcome through GitHub Issues.
Please include a minimal reproducible example, `sessionInfo()`, the selected
reference version, and checksums recorded in the ProReT run manifest.

Contributions should be submitted as pull requests. New behavior requires
documentation and tests. Changes to numerical output must describe the
scientific rationale, preserve backward compatibility where possible, and add
or update a frozen toy regression test. Run `devtools::test()` and
`rcmdcheck::rcmdcheck()` before submission.

By contributing, you agree that your contribution is released under the MIT
license used by this repository.
