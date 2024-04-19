```@meta
CurrentModule = FerriteInterfaceElements
```

# FerriteInterfaceElements.jl

Welcome to the documentation for [FerriteInterfaceElements](https://github.com/Ferrite-FEM/FerriteInterfaceElements.jl)! FerriteInterfaceElements adds functionality to the
finite element toolbox [Ferrite](https://github.com/Ferrite-FEM/Ferrite.jl) for working
with interface elements (also known as cohesive elements).

## How the documentation is organized

The documentation assumes that you are already familiar with the basic usage of Ferrite.
If not, you should first take a look at the [Ferrite documentation](https://ferrite-fem.github.io/Ferrite.jl/dev/).
Here, only the additional tools are explained.

After a basic introduction on this side, the document is organized as follows:

 - [**Tutorials**](tutorials/index.md) are documented examples which guide you
   through the process of solving partial differential equations with interface terms.
 - [**Reference**](reference/index.md) contains the technical API reference of functions and
   methods (e.g. the documentation strings).

## Getting started

As a new user of FerriteInterfaceElements it is suggested to read the introduction on this side
and then start working with the tutorials before using FerriteInterfaceElements to tackle the 
specific equation you ultimately want to solve.

### Getting help

If you have questions about FerriteInterfaceElements it is suggested to use the `#ferrite-fem` channel on the
[Julia Slack](https://julialang.org/slack/), or the `#Ferrite.jl` stream on
[Zulip](https://julialang.zulipchat.com/). Alternatively you can use the [discussion
forum](https://github.com/Ferrite-FEM/FerriteInterfaceElements.jl/discussions) on the GitHub repository.

### Installation

To use FerriteInterfaceElements you first need to install Julia, see <https://julialang.org/> for details.
Installing FerriteInterfaceElements can then be done from the Pkg REPL; press `]` at the `julia>` promp to
enter `pkg>` mode:

```
pkg> add FerriteInterfaceElements
```

This will install FerriteInterfaceElements and all necessary dependencies. Press backspace to get back to the
`julia>` prompt. (See the [documentation for Pkg](https://pkgdocs.julialang.org/), Julia's
package manager, for more help regarding package installation and project management.)

Note, that you also need to install Ferrite, which can be done in the same way.

Finally, to load Ferrite and FerriteInterfaceElements, use

```julia
using Ferrite, FerriteInterfaceElements
```

You are now all set to start using FerriteInterfaceElements!

### Introduction to interface elements

Interface elements can be used to allow for jumps in fields on interfaces between
finite elements. They can be considered as elements with zero thickness and
two faces (one for each side of the interface).

In FerriteInterfaceElements, their implementation is based on the idea to combine
two embedded elements to represent the to sides of an interface. At the current state,
it is only supported to use the same interpolation on both sides. To distinguish
the two sides, they are referred to as *here* and *there*. The jump of a
field is then defined as from *here* to *there*: *field there - field here*.

Another key definition is the reference shape of interface elements. In
FerriteInterfaceElements, the shape is chosen to be of the same dimension as the
surrounding bulk elements. The specific shape is then determined from the shape
of the underlying base shapes:

| reference shape of faces | shape of the interface element |
|:------------------------ |:------------------------------ |
| `RefLine`                | `RefQuadrilateral`             |
| `RefTriangle`            | `RefPrism`                     |
| `RefQuadrilateral`       | `RefHexahedron`                |

These definitions have been made to use as much of the existing functionalities in
Ferrite as possible on the types for working with interface elements:

* `InterfaceCell`
* `InterfaceCellInterpolation`
* `InterfaceCellValues`

These types can be used similar to the corresponding types in Ferrite. However, they
come with new functions for evaluations:

* `getdetJdV_average`
* `shape_value_jump`
* `shape_value_average`
* `function_value_jump`
* ...