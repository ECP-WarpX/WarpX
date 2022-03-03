.. _development-faq:

FAQ
===

This section lists frequently asked developer questions.


What is ``0.0_rt``?
-------------------

It's a C++ `floating-point literal <https://en.cppreference.com/w/cpp/language/floating_literal>`__ for zero of type ``amrex::Real``.

We use `literals <https://en.cppreference.com/w/cpp/language/expressions#Literals>`__ to define constants with a specific type, in that case the zero-value.
There is also ``0.0_prt``, which is a literal zero of type ``amrex::ParticleReal``.
In std C++, you know: ``0.0`` (literal ``double``), ``0.0f`` (literal ``float``) and ``0.0L`` (literal ``long double``).
We do not use use those, so that we can configure floating point precision at compile time and use different precision for fields (``amrex::Real``) and particles (``amrex::ParticleReal``).

You can also write things like ``42.0_prt`` if you like to have another value than zero.

We use these `C++ user literals <https://en.cppreference.com/w/cpp/language/user_literal>`__ (`[1] <https://github.com/AMReX-Codes/amrex/pull/577>`__, `[2] <https://github.com/AMReX-Codes/amrex/pull/578>`__, `[3] <https://github.com/AMReX-Codes/amrex/pull/869>`__), because we want to avoid that double operations, i.e., ``3. / 4.``, implicit casts, or even worse integer operations, i.e., ``3 / 4``, sneak into the code base and make results wrong or slower.


Do you worry about using ``size_t`` vs. ``uint`` vs. ``int`` for indexing things?
---------------------------------------------------------------------------------

`std::size_t <https://en.cppreference.com/w/cpp/types/size_t>`__ is the C++ unsigned int type for all container sizes.

Close to but not necessarily ``uint``, `depends on the platform <https://en.cppreference.com/w/cpp/language/types>`__.
For "hot" inner loops, you want to use ``int`` instead of an unsigned integer type. Why? Because ``int`` has no handling for overflows (it is intentional, undefined behavior in C++), which allows compilers to vectorize easier, because they don't need to check for an overflow every time one reaches the control/condition section of the loop.

C++20 will also add support for `ssize <https://en.cppreference.com/w/cpp/iterator/size>`__ (signed size), but we currently require C++17 for builds.
Thus, sometimes you need to ``static_cast<int>(...)``.


What does ``std::make_unique`` do?
----------------------------------

`make_unique <https://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique>`__ is a `C++ factory method <https://refactoring.guru/design-patterns/factory-method/cpp/example>`__ that creates a ``std::unique_ptr<T>``.

Follow-up: ``Why use this over just *my_ptr = new <class>``?

Because so-called `smart-pointers <https://en.cppreference.com/book/intro/smart_pointers>`__, such as ``std::unique_ptr<T>``, do delete themselves automatically when they run out of scope.
That means: no memory leaks, because you cannot forget to ``delete`` them again.


Why name header files ``.H`` instead of ``.h``?
-----------------------------------------------

This is just a :ref:`convention that we follow through the code base <developers-contributing-style-conventions>`, which slightly simplifies what we need to parse in our various build systems.
We inherited that from AMReX.
Generally speaking, C++ file endings can be arbitrary, we just keep them consistent to avoid confusion in the code base.

To be explicit and avoid confusion (with C/ObjC), we might change them all to ``.hpp`` and ``.cpp``/``.cxx`` at some point, but for now ``.H`` and ``.cpp`` is what we do (as in AMReX).


What does const int ``/*i_buffer*/`` mean in argument list?
-----------------------------------------------------------

This is often seen in a derived class, overwriting an interface method.
It means we do not name the parameter because we do not use it when we overwrite the interface.
But we add the name as a comment ``/* ... */`` so that we know what we ignored when looking at the definition of the overwritten method.
