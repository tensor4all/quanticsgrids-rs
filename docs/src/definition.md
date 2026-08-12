# Definition

We first introduce a base-`B` representation with `B = 2, 3, 4, ...`.
Throughout this guide, quantics digits and grid indices are 0-based
(0-indexed), matching the Rust convention. QuanticsGrids.jl is 1-based;
when porting Julia scripts, subtract 1 from grid indices and quantics
digits at the call boundary.

We represent a non-negative integer `X >= 0` as

$$
X = \sum_{i=1}^{R} x_i B^{R-i},
$$

where each digit `x_i` satisfies `0 <= x_i < B` and `R` is the number of
digits. In this crate, the base-`B` representation of `X` is stored as the
vector

$$
[x_1, \ldots, x_R].
$$

For multiple variables, the crate supports fused and interleaved unfolding
schemes. For three variables `X`, `Y`, and `Z`, suppose their base-`B`
representations are

$$
[x_1, \ldots, x_R], \quad [y_1, \ldots, y_R], \quad [z_1, \ldots, z_R].
$$

The interleaved representation is

$$
[x_1, y_1, z_1, x_2, y_2, z_2, \ldots, x_R, y_R, z_R].
$$

The fused representation is

$$
[\alpha_1, \alpha_2, \ldots, \alpha_R],
$$

where

$$
\alpha_i = x_i + B y_i + B^2 z_i
$$

and therefore

$$
0 \le \alpha_i < B^3.
$$

In fused ordering, the `x` digit runs fastest at each digit level. This
generalizes to any number of variables.
