Oil Price Shocks in a Small Open Economy New
Keynesian Model with Input-Output Linkages∗

Agustin Diaz † Luigi Durand Jorge Fornero

Benjam´ın Garc´ıa Mario Giarda Josefina Henriquez

June 8, 2026

Abstract

We build a New Keynesian small open economy model (nk-iosoe) with input-output
production networks calibrated to Chile’s 12-sector economy. Sectors combine labor,
domestic intermediate inputs determined by the empirical Chilean IO matrix, and im-
ported goods. Price stickiness is heterogeneous across sectors and labor reallocation
involves quadratic adjustment costs. We estimate structural parameters by the Simu-
lated Method of Moments, matching sectoral output, price, and employment volatilities
in Chilean data from 2003 to 2023. The production network is a quantitatively impor-
tant amplifier: the IO multiplier raises the inflationary impact of a 10% world oil price
increase by a factor of 2.5 relative to the direct cost channel, with network propaga-
tion accounting for 60 percent of the total marginal cost increase. Heterogeneous price
stickiness interacts with the IO structure to generate cross-sectoral dispersion in price
responses that representative-sector models cannot replicate.
JEL codes: E31, E52, F41, L16.
Keywords: Input-output networks, small open economy, New Keynesian model, price
rigidity, oil price shocks, Chile.

∗Aknowledges:..
†Central Bank of Chile, adiaz@bcentral.cl

1

Introduction

Understanding how shocks propagate through an economy has long been at the center of

macroeconomics, but the role of production networks in shaping that propagation has only

recently received rigorous attention.

Input-output linkages create two-way relationships

between sectors: cost shocks in upstream suppliers raise marginal costs in downstream buy-

ers, while demand disruptions in large customers reduce output in their suppliers.

In an

economy with nominal rigidities, these upstream and downstream linkages interact with

sectoral price stickiness to generate rich patterns of inflation and output dynamics that a

representative-sector model cannot capture. This interaction is particularly consequential in

small open economies, where traded commodities—most importantly oil—enter production

as intermediate inputs at highly heterogeneous rates across sectors, and where exchange-rate

movements propagate through the same production network.

This paper studies these questions for Chile. Chile is a natural laboratory:

it is a

small open economy with well-documented sectoral heterogeneity in price stickiness (Pasten

et al., 2021), a rich input-output structure, and significant exposure to commodity price

shocks—particularly oil, which constitutes a large share of intermediate imports in Trans-

port, Mining, and Manufacturing. We build a New Keynesian small open economy model

with N = 12 sectors whose linkages are fully pinned down by the Chilean 2021 input-output

table, calibrate it to Chilean data from 2003 to 2023, and use it to study the transmission

of preference, technology, monetary policy, and oil price shocks.

What we do. Our model, which we term nk-iosoe (New Keynesian Input-Output Small

Open Economy), extends Romero (2023) and Ferrante et al. (2023) to a setting with (i)

a fully estimated 12-sector IO structure for Chile, (ii) sector-specific Rotemberg price ad-

justment costs calibrated from micro price data, (iii) a disaggregated oil-versus-non-oil split

of intermediate imports at the sector level, and (iv) a three-layer consumption aggregator

distinguishing goods from services, sectoral varieties, and home from foreign production. We

estimate the remaining structural parameters—production elasticities, shock processes, and

the export demand elasticity—by Simulated Method of Moments, targeting 15 unconditional

second moments of Chilean quarterly data: the cross-sectional dispersion and rank ordering

1

of sectoral output, price, and employment volatilities, along with aggregate GDP volatility,

inflation volatility, their correlation, the real exchange rate volatility, and the average trade

balance ratio.

Main findings. Three results stand out. First, the production network is a quantita-

tively important amplifier of shocks. Taking the oil price shock as our leading example,

the Leontief-inverse multiplier amplifies the aggregate inflationary impact of a 10% world

oil price increase by a factor of 2.5 relative to the direct cost channel alone—IO network

propagation accounts for 60 percent of the total marginal cost increase. Sectors with mod-

erate own-oil exposure but high upstream dependence on oil-intensive suppliers—such as

Agriculture and Business Services—experience large second-round cost increases that the di-

rect channel alone would miss. Second, heterogeneous price stickiness interacts non-trivially

with the IO structure. Sectors with stickier prices absorb a larger fraction of cost shocks

as markup compression in the short run, dampening their own price response while am-

plifying the cost push transmitted to downstream buyers. This interaction reshuffles the

cross-sectoral ranking of inflationary impacts relative to what either the network structure

or the price-stickiness distribution would predict independently. Third, the SOE dimension

matters for monetary policy. Exchange-rate depreciation in response to a negative demand

shock amplifies the import cost channel, while the Taylor rule faces a trade-off between sta-

bilizing the output contraction and controlling the imported inflation component, a trade-off

that is absent in closed-economy IO models.

Related literature. This paper contributes to several strands of the literature. The

theoretical underpinning of production networks in macroeconomics was developed by ?

and ?, who show that the Leontief inverse governs first-order aggregate effects of sectoral

shocks. Ferrante et al. (2023) introduces IO linkages into a New Keynesian model and shows

they amplify inflationary supply shocks; our model generalizes their framework to an open

economy with sector-specific price stickiness. Romero (2023) studies exchange-rate pass-

through in an SOE IO model and shows that the network structure matters for how nominal

depreciations transmit to domestic inflation across sectors; we build on his specification and

2

extend it with a full calibration and shock analysis for Chile. On the IO-commodity price

nexus, Romero (2022) and Silva et al. (2024) show that up- and downstream positions in the

IO network determine a sector’s vulnerability to commodity shocks. Rubbo (2023) derives

the sectoral Phillips curves implied by IO networks and shows that monetary policy faces

a richer trade-off in multi-sector economies. On the Chilean economy specifically, Pasten

et al. (2021) provide the micro-data evidence on sectoral price stickiness that we use for

calibration. Our contribution is to integrate these building blocks into a single calibrated

framework and use it for a quantitative analysis of shock transmission in Chile.

Paper structure. Section 2 presents the model. Section 3 describes the calibration and

data. Section 4 analyzes the transmission of four structural shocks, with particular emphasis

on the oil price shock and the decomposition of direct-cost versus network-amplification

effects. Section 5 concludes. Technical derivations and additional tables are in the Appendix.

2 The Model

The economy consists of N = 12 sectors in a small open economy. Sectors produce using

three inputs—labor, domestic intermediate inputs, and imported goods—whose combination

is described by a CES production function. The intermediate input bundle is itself a CES ag-

gregate of outputs from all N sectors, with weights determined by the empirical input-output

matrix. Sectors are subject to heterogeneous Rotemberg price adjustment costs calibrated

from micro price data, and labor reallocation across sectors involves quadratic adjustment

costs. Production can be directed to household consumption, export, or intermediate use in

other sectors. The representative household purchases goods from all sectors organized in a

three-layer CES structure distinguishing goods from services, individual sectors within each

category, and home from foreign varieties, following Gal´ı and Monacelli (2005). Households

save in risk-free local and foreign bonds subject to a debt-elastic foreign interest rate that

generates a well-behaved external debt position.

3

2.1 Firms

Sector i final good producer.

In each sector i, a representative competitive producer

aggregates the output of a continuum of intermediate monopolistically competitive firms

index by s:

Y i
t =





Z 1

0



ϵ
ϵ−1

Y i
t (s)

ϵ−1
ϵ ds



where ϵ is the elasticity of substitution between varieties within a sector. Then, the solution

to the competitive producer’s problem implies the following demand curve for differentiated

products in each sector:

Y i
t (s) =

  P i
t (s)
P i
t

!−ϵ

Y i
t

Sector i intermediate goods producers.

In each sector, a continuum of firms supply

differentiated products to the representative competitive producer subject to price adjust-

ment costs. The differentiated products are produced according to the following production

function,

t (s) = Ai
Y i
t

1
mi M i
ϵY

t (s)

α

ϵY −1
ϵY + α

1
vi V i
ϵY

t (s)

ϵY −1
ϵY + (1 − αmi − αvi)

1

ϵY Li

t(s)

ϵY −1
ϵY

! ϵY

ϵY −1

,

where ϵY is the elasticity of substitution between production inputs, αmi determines the

materials share, and αvi imported goods share. The intermediate inputs are a bundle of the

outputs of the N sectors of the economy:

M i

t (s) =





N
X

j=1

1
ϵm

i,j (M i
Γ

j,t(s))



ϵm
ϵm−1



,

ϵm−1
ϵm

where ϵm is the elasticity of substitution between the different sectors intermediate inputs.
The input-output matrix parameter Γi,j (where PN

j=1 Γi,j = 1, ∀ i) determines the share of

inputs that sector i purchases from sector j. Note that this uses the purchasing sector as

the first index, which is the reverse of the standard IO table convention shown in Table 8.

4

Input demands. The cost—intratemporal minimization problem of firms is given by

{M i

jt(s)}N

t (s),Li

t(s)

min
j=1,V i

N
X

j=1

jt M i
P H

jt(s) + P v

itV i

t (s) + W i

t Li

t(s)

s.t.

t (s) =Ai
Y i
t

1
mi M i
ϵY
α

t (s)

ϵY −1
ϵY + α

1
vi V i
ϵY

t (s)

ϵY −1
ϵY + (1 − αmi − αvi)

1

ϵY Li

t(s)

! ϵY

ϵY −1

ϵY −1
ϵY



M i

t (s) =



1
ϵm

i,j (M i
Γ

j,t(s))



ϵm
ϵm−1



ϵm−1
ϵm

N
X

j=1

The FOCs are

∀ j = 1, . . . , N,

i = 1, . . . , N

(1)

P H
jt
Pt

P M i
t
Pt

P iV
t
Pt

P Li
t
Pt

=M C i

t(s)

αmi

=M C i

t(s)

αmi

! 1
ϵY

! 1
ϵY

Y i
t (s)
M i
t (s)

Y i
t (s)
M i
t (s)

=M C i

t(s)

αvi

! 1
ϵY

Y i
t (s)
V i
t (s)

! 1
ϵm

Γij

M i
M i

t (s)
jt(s)

∀ i = 1, . . . , N

∀ i = 1, . . . , N

=M C i

t(s)

(1 − αmi − αvi)

! 1
ϵY

Y i
t (s)
Li
t(s)

∀ i = 1, . . . , N

(2)

(3)

(4)

Notice that the price of materials is given by the home price of the material, P H
jt .

Given the CES aggregator, the cost minimization problem implies the following price

indices for intermediates inputs in sector i, which comprises only home (H) goods:

P M i

t =





N
X

j=1

1
ϵm

i,j (P H
Γ

jt )1−ϵm



1
1−ϵm



Given the price index for intermediate inputs, the price of imported inputs, and wages

in sector i, the marginal cost of production in sector i is:

M C i

t =


αmi(P M i

t

1
Ai
t

)1−ϵY + αvi(P v

t )1−ϵY + (1 − αmi − αvi)(P Li

t )1−ϵY



1
1−ϵY



Log-linearizing the marginal cost expression and the material price index around the

5

symmetric steady state (where all prices equal one and cost-share weights coincide with the

calibrated factor shares) yields

dmci

t = αmi

N
X

Γij bpH

jt

j=1
{z
network (via input prices)

|

}

+ αvi αOil
i
{z
direct oil

|

bpO
t
}

+ αvi(1 − αOil

|

{z
non-oil imports

i ) bpv
t
}

+ (1 − αmi − αvi) bwt
{z
}
labor

|

− ˆAi
t,

(5)

where bxt ≡ log xt − log ¯x denotes the log-deviation from steady state. Equation (5) decom-
poses marginal cost changes into four economically distinct channels. The network term

shows that sector i’s costs rise when any upstream supplier j raises its price P H

j , weighted

by the IO share Γij; the strength of this channel is governed by the materials share αmi. The

direct oil term isolates the first-round cost push from the domestic oil price P O

t , scaled by

the imported-input share αvi and the sector-specific oil intensity αOil

i

. The non-oil import

term captures the exchange-rate-driven component of non-oil intermediate import costs. The

labor term enters with weight (1 − αmi − αvi), the labor share in production costs. Because

all four components use endogenous GE prices, equation (5) provides an exact accounting

of marginal cost changes in the linearized model; it is used in Section 4.3 to decompose the

sectoral transmission of oil price shocks.

We define the real exchange rate as

Qt ≡

EtP ∗
t
Pt

,

(6)

where P ∗
t

is the foreign CPI and Pt is the domestic CPI. An increase in Qt represents a real

depreciation. Under our baseline assumption of complete exchange-rate pass-through, the

domestic-currency price of imported goods moves one-for-one with the nominal exchange

rate, so that P V

t = QtP V ∗

t

.

6

Intermediate Producers: Price setting. Given the marginal cost just derived, home

firms set prices subject to non-pecuniary, quadratic cost adjustment. The firm problem is:

  P H
it (s)
Pt

!−ϵ Y i
t
Pt

(P H

it (s) − M C i

t) −

κi
2

  P H
it (s)
P H

it−1

!2

− 1

P H
it

Y i
t
Pt

t (P H
V i

it−1(s)) = max
P H
it (s)
"

+ Et

Mt+1V i

t+1(P H

#
,
it (s))

where κi is the sector specific adjustment cost, Mt+1 is the stochastic discount factor of

the representative household. The solution to the previous problem is the sector-level New

Keynesian Phillips Curve:

1 − ϵ + ϵ

M C i
t
P H
it

− κi(ΠH

it − 1)ΠH

it + κiEt

Mt+1ΠH

it+1(ΠH

it+1 − 1)

!

Y i
t+1
Y i
t

= 0

where Πt ≡ Pt/Pt−1 is the gross aggregate CPI inflation rate and ΠH
it

is the gross nominal

sectoral inflation rate, defined as the quarter-over-quarter ratio of the nominal home price
ˆP H
it :

ΠH

it ≡

ˆP H
it
ˆP H

it−1

.

(7)

In the model’s implementation it is convenient to express prices in units of the aggregate

CPI. We therefore define the normalized home price

P H

it ≡

ˆP H
it
Pt

,

(8)

so that P H

it measures the home price of sector i relative to the aggregate price level. In terms

of this normalized price, nominal sectoral inflation decomposes as

ΠH

it =

ˆP H
it
ˆP H

it−1

=

P H
it
P H

it−1

Πt,

(9)

linking nominal sectoral inflation to the change in the relative home price and the common

aggregate inflation rate. The decomposition in (9) is consequential for impulse response

analysis: a model-based IRF for P H

it traces only the movement in the relative home price.

Nominal sectoral inflation ΠH

it —the object appearing in the NKPC and governing price-

7

setting incentives—is recovered by adding back aggregate CPI inflation Πt.

2.2 Importable Goods and Foreign Demand for Produced Goods

Our small open economy setting assumes that some produced goods are exported, sectors

use imported goods as production inputs, and consumers purchase imported goods directly.

Exportable good. Sectors i = 1, ..., N are demanded from abroad by a competitive firm

which combines them as inputs to produce an exportable good. The technology of this sector

is given by

Y x
t =δx

  N
X

i=1

1
it (Y x
ϵX
it )
ψ

ϵX −1
ϵX

!

ϵX
ϵX −1

(10)

where δx is a constant term, PN

j=1 ψj = 1 holds, and ψit controls the share of output exported

by sector i and can be a foreign demand shock for sector i production. The price of exports

and the demand of goods from each sector are given by

P x

t =





N
X

j=2

ψitP 1−ϵX
jt



1
1−ϵX



and Y x

jt = ψj

!ϵX

  P x
t
Pjt

Y x
t .

Additionally, the foreign demand for the exportable good takes the form Y x

t = ωx (cid:16) P x

t
EtP ∗
t

(11)

(cid:17)−η∗

Y ∗
t ,

where P ∗

t and Y ∗

t denote the foreign CPI and the foreign level of output, respectively.

Imported inputs: oil and non-oil. Production imports V i
t

in each sector consist of an

oil component V iO

t

and a non-oil component V iN

t

, bundled via a CES aggregator:

(cid:18)(cid:16)

V i
t =

αiV
N

(cid:17) 1

ϵV V iN
t

ϵV −1
ϵV +

(cid:16)

αiV
O

(cid:17) 1

ϵV V iO
t

ϵV −1
ϵV

(cid:19)

ϵV
ϵV −1

,

(12)

where αiV

O is the oil share in sector i’s import bundle, αiV

N = 1 − αiV

O is the non-oil share, and

ϵV is the elasticity of substitution between oil and non-oil imports. The sector-level oil shares

αiV

O are taken from the 2021 Chilean Input-Output tables (Table 7). Cost minimization yields

8

the input demands

t = αiV
V iO
O

t = αiV
V iN
N

!ϵV

!ϵV

  P iV
t
P O
t
  P iV
t
P V
t

V i
t ,

V i
t ,

and the sector-specific import price index

(cid:16)

P iV

t =

O P O 1−ϵV
αiV
t

+ αiV

N P V 1−ϵV
t

(cid:17) 1

1−ϵV .

(13)

(14)

(15)

The non-oil import price is P V

t = QtP V ∗

t

as described above. The oil price is P O

t = QtP O∗

t

,

where the international oil price P O∗

t

follows an AR(1) process:

log

!

  P O∗
t
P O∗

= ρpo log

!

 P O∗
t−1
P O∗

+ σpo εo
t .

The trade balance therefore separates into oil and non-oil components:

T Bt = P X

t Xt − P V
t

(cid:16)P

i V iN

t + C F
t

(cid:17)

X

− P O
t

V iO
t

.

i

(16)

(17)

Imported consumption goods. We assume there is an homogeneous non-oil imported

good demanded directly for household consumption, denoted by C iF
t

. Its price is P V
t

(the

non-oil import price defined above). Household demand for imported variety i is given by

the Armington structure in equations (??).

Our baseline model assumes complete exchange-rate pass-through to import prices. The

domestic price of non-oil imports is P V

t = QtP V ∗

t

, where P V ∗

t

is the world price of the

importable good (expressed in foreign currency), which follows the AR(1) process

log

!

  P V ∗
t
P V ∗

= ρpv log

!

  P V ∗
t−1
P V ∗

+ σpv εpv
t ,

(18)

with estimated parameters ρpv and σpv (Table 3). Under this specification, changes in the

nominal exchange rate pass through one-for-one into domestic import prices in each period.

9

2.3 Households

The household dynamic problem is:

Vt(Bt, B∗

t ) = max
Ct,Nt

C 1−γ
t
1 − γ

− χt

N 1+ψ
t
1 + ψ

+ βVt+1(Bt+1, B∗

t+1)

s.t PtCt + Bt+1 + EtB∗

t+1 = WtNt + RtBt + EtR∗

t B∗

t + divt,

where divt are the profit of the monopolistically competitive firms, labor agencies, and profits

of the importer sector. Bt and B∗

t are local and foreign nominal bond holdings, χt is the

disutility of labor supply that is allowed to vary over time, and WtNt represents households’

labor income. Then the household first order conditions are:

ξtC −γ

t = βEt


ξt+1C −γ

t+1Rt+1


,

1
Πt+1

C −γ
t

Wt
Pt

ξtC −γ

t = βEt

= χtN ψ
t ,


ξt+1C −γ

t+1R∗

t+1


,

ΠE
t+1
Πt+1

where the stochastic discount factor is Mt+1 = β

!−γ

Ct+1
Ct

The preference shock ξt enters the Euler equation multiplicatively and follows the AR(1)

process:

ξt = (1 − ρξ) + ρξξt−1 + σξεξ
t ,

where εξ

t ∼ N (0, 1). A positive shock ξt > 1 raises the current marginal utility of consump-

tion, generating a demand expansion with positive output-inflation comovement. Parameters

ρξ and σξ are estimated jointly with the other shock processes (Table 3).

Consumption aggregators.

In this multisector economy, total consumption Ct is a bun-

dle of several aggregators: goods-services, the N sectors, and home-foreign produced goods.

First, the goods and services aggregator with g = {S, G} combines consumption of ser-

10

vices and goods in a Cobb-Douglas form,

Ct =

 C G
t
ΩG

!ΩG   C S
t
ΩS

!ΩS

,

with ΩG + ΩS = 1. This Cobb-Douglas specification implies that households allocate fixed

expenditure shares ΩG and ΩS to goods and services, respectively, independent of relative

prices.

Second, within each g, consumers purchase all sectors i and combine them with a Cobb-

Douglas aggregator

C g

t =

N
Y

i=1

(C g

it)ωg

i

g = {S, G},

with g-specific weights ωg

i > 0 satisfying PN

i=1 ωg

i = 1. This Cobb-Douglas specification

implies constant within-category expenditure shares.

And third, each C g

it is provided by a home H or a foreign F producer. Following Gal´ı

and Monacelli (2005), these are bundled with

C g

it =

(cid:18)
ϱiC Hg
it

σH −1
σH + (1 − ϱi)C F g
it

σH −1
σH

(cid:19) σH
σH −1

i = {1, ..., N },

g = {S, G}

(19)

where ϱi < 1 and represents the home bias of consumption in good (i, g).

Consumption allocation problem and FOCs. The first order conditions of the first

layer are given by

with price index given by

C g

t = Ωg

Pt
P g
t

Ct,

g = {S, G},

Pt =

(cid:16)

P G
t

(cid:17)ΩG (cid:16)

(cid:17)ΩS .

P S
t

11

The demands for each good within category g = {S, G} are

C g

it = ωg

i

P g
t
Pit

C g
t

i = {1, ..., N }, g = {S, G}.

The price index of goods and services is

P g

t =

P ωg
it .

i

N
Y

i=1

Then, among each (i, g), we solve for the home and foreign origin of goods as

C Hg

it = ϱσH

i

 P H
it
P g
it

!−σH

C g
it,

i = {1, ..., N }, g = {S, G},

C F g

it = (1 − ϱi)σH

  P F
it
P g
it

!−σH

C g
it,

i = {1, ..., N }, g = {S, G}.

The (average) price of good i by the side of households is given by

(cid:16)

P g

it =

ϱσH
i P H
it

1−σH + (1 − ϱi)σH P F
it

1−σH

(cid:17) 1

1−σH

Notice that P G

it = P S

it because ϱi is assumed to be the same for all g’s. Thus, we can

confidently denote these prices as Pit. Additionally, we assume all imported goods are

homogeneous, which implies P F

it = P v
t

∀ i = 1, ..., N .

2.4 Labor Markets

In each sector, labor is supplied to the monopolistically competitive firm by a representative

labor agency that hires labor from the representative household. Agencies face a quadratic

adjustment cost to switching labor. Then, the agency’s problem is:

W i

t (Li

t−1) = max
t,גi
Li
jt

P L,i
t
Pt

Li

t −


1 +

Li
t

Wt
Pt

c
2

  Li
t
Li

t−1

!2

"

− 1

 + Et

Mt+1W i

#
t+1(Li
t)

where c is the hiring cost. From this problem, we can recover the sectoral labor demand:

12

P L,i
t
Pt

=

Wt
Pt

+

"

Wt
Pt

  c
2

  Li
t
Li

t−1

!2

− 1

+ c

−Et

Mt+1

Wt+1
Pt+1

  Li
t+1
Li
t

c

− 1

− 1

! Li
t
Li

t−1





!2#

t−1

  Li
t
Li
! Li
t+1
Li
t

current and future expected labor costs introduce a wedge between the aggregate wage

and the price of labor in each sector. The wedge generates dividends that are distributed to

the household.

2.5 Monetary Policy and Market Clearing

Market clearing in each sector requires that domestic supply equals the sum of domestic

consumption, exports, and intermediate demand from all other sectors:

Yjt = C H

jt + Xjt +

N
X

i=1

Mijt.

(20)

The trade balance in nominal terms combines exports and the full import bill, separating oil

from non-oil imports:

T Bt = P X

t Xt − P V
t

(cid:16)P

i V iN

t + C F
t

(cid:17)

X

− P O
t

V iO
t

.

i

Gross domestic product is the sum of consumption and the trade balance:

GDPt = Ct + T Bt.

(21)

(22)

Local bonds are in zero net supply (Bt = 0). The external position evolves according to

EtB∗

t = T Bt + EtR∗

t−1B∗

t−1.

(23)

13

The domestic nominal interest rate is set by the central bank via a Taylor rule with smooth-

ing:

Rt = (1 − ρi)

1
β

+ ρiRt−1 + (1 − ρi)ϕπ(Πt − 1) + εi
t,

(24)

where ρi = 0.74 is the interest-rate-smoothing coefficient, ϕπ = 2.5 governs the response to

inflation deviations from target, and εi

t is a monetary policy shock.

To close the model, the foreign interest rate includes a debt-elastic risk premium:

R∗

t = Rw

(cid:20)
ξb
t × exp

(cid:18)
¯b −

EtB∗
t
GDPt

(cid:19)(cid:21)

,

(25)

where ξb > 0, Rw
t

is the world risk-free rate, and ¯b is the steady-state net foreign asset

position (calibrated to match Chile’s average trade balance). This specification ensures a

well-defined, stationary external debt position: when foreign debt rises above its long-run

target, the cost of borrowing rises, discouraging further accumulation.

The aggregate labor market clears when total effective labor supply (inclusive of adjust-

ment costs) equals household labor Nt:

N
X

i=1

Li
t

1 +

c
2

  Li
t
Li

t−1

!2!

− 1

= Nt.

3 Calibration and Data

This section describes the calibration strategy and data sources for the model parameters. We

organize parameters into three categories based on their source: (i) parameters set from the

literature, (ii) parameters constructed directly from data, and (iii) parameters that require

calibration or estimation.

14

Table 1. Model Calibration: Parameters from Literature

Parameter Description

Value

Source

Production Technology

ϵ

Elast. subst., varieties

10

Ferrante et al. (2023)

Price Setting

κi

Price adj. cost, sector i

Table 4

Romero (2023)

Household

β

γ

ψ

χ

Discount factor (quarterly)

0.986

Standard

Inv. intertemporal elast.

Inv. Frisch elasticity

Labor disutility weight

2

1

1

Ferrante et al. (2023)

Ferrante et al. (2023)

Ferrante et al. (2023)

Monetary Policy

ϕπ

ρi

Taylor rule coef. (inflation)

2.5

Literature

Interest rate smoothing

0.74

Literature

Table 2. Model Calibration: Parameters from Data

Parameter Description

Value

Source

Production Technology

Continued on next page

15

Table 2 – Continued from previous page

Parameter Description

Value

Source

αmi

αvi

Γi,j

Consumption

ΩG

ΩS

ωg
i

ωs
i

ϱi

Material share, sector i

Table 4

I-O Tables (Subsection 3.5.2)

Import share, sector i

Table 4

I-O Tables (Subsection 3.5.2)

IO matrix (model: purchaser × supplier)

Table 8

I-O Tables (Subsection 3.5.1)

Goods weight in consumption

0.57

National Accounts (Subsection 3.5.3)

Services weight in consumption

0.43

National Accounts (Subsection 3.5.3)

Goods sector i weight

Table 5

National Accounts (Subsection 3.5.3)

Services sector i weight

Table 5

National Accounts (Subsection 3.5.3)

Home bias, sector i

Table 4

National Accounts (Subsection 3.6)

Small Open Economy

P ∗
t

Y ∗
t

Π∗
t

Rw

χX
i

Foreign price level (SS)

Foreign output (SS)

Foreign inflation (SS, quarterly)

World interest rate (SS, quarterly)

1

1

1.00

1.045

Normalized

Normalized

Normalized

Data

Export share, sector i

Sector-specific

I-O Tables

16

Table 3. Model Calibration: Calibrated and Estimated Parameters

Parameter Description

Value

Method

Production Technology

Elast. subst., production inputs

1.484

Elast. subst., materials (IO sectors)

0.051

SMM

SMM

εY

εm

Price Setting

κv

Import price adjustment cost

1,382

SMM

Consumption

σH

Labor

c

ψl

Elast. subst., home-foreign goods

0.999

Calibrated

Aggregate labor reallocation cost

0.025

SMM

Sectoral labor adj. cost

0

Calibrated

Small Open Economy

η∗

ωx

Financial

ξb

¯b

Export demand elasticity

Foreign demand scale

0.500

1

SMM

Calibrated

Debt elasticity (risk premium)

0.001

Calibrated

Steady-state debt-to-GDP

−0.57 Endogenous (from TB target)

Continued on next page

17

Table 3 – Continued from previous page

Parameter Description

Value

Method

Shocks: Goods-Services Preference

ρω

σω

Preference shock, AR(1) coef.

Preference shock, std. dev.

0.591

0.138

Shocks: Technology (TFP)

ρA

¯σA

TFP shock, AR(1) coef. (common)

0.446

TFP shock, std. dev. (cross-sector mean)

0.054

Shocks: Monetary Policy

SMM

SMM

SMM

SMM

σi

Monetary shock, std. dev.

0.001

Calibrated

Shocks: World Import Price

ρP ∗

σP ∗

Non-oil import price, persistence

Non-oil import price, std. dev.

0.515

0.056

Shocks: Foreign Demand

ρξ

σξ

Foreign demand shock, persistence

0.462

Foreign demand shock, std. dev.

0.042

SMM

SMM

SMM

SMM

Literature-based parameters (Table 1) Standard macroeconomic parameters follow

Ferrante et al. (2023): the elasticity of substitution between varieties (ϵ = 10), risk aversion

coefficient (γ = 2), inverse Frisch elasticity (ψ = 1), and labor disutility weight (χ = 1).

The discount factor (β = 0.986) is calibrated to match Chile’s average long-run real interest

18

rate. Monetary policy rule coefficients (ϕπ = 2.5, ρi = 0.74) are consistent with estimates

for small open economies in the literature.

Nominal rigidities. We discipline sectoral price stickiness using the Calvo non-adjustment

probabilities from Romero (2023). Specifically, we take two values—one for goods and one

for services—and assign θg to each goods-producing sector and θs to each services-producing

sector (with θj denoting the probability that a firm in sector j does not reoptimize its price in

a given period). For computational convenience, we implement Rotemberg (quadratic) price

adjustment costs. To preserve the degree of nominal rigidity implied by the Calvo estimates,

we convert θj into a sector-specific Rotemberg parameter κj by matching the slope of the

log-linear New Keynesian Phillips curve around a zero-inflation steady state:

κj = (ε − 1)

θj
(1 − θj) (1 − βθj)

,

where ε is the within-sector elasticity of substitution and β is the discount factor. This

mapping ensures that each sector’s inflation dynamics under Rotemberg pricing replicate, to

first order, the price stickiness implied by the Calvo probabilities in Romero (2023).

Data-based parameters (Table 2) Production structure parameters are derived from

Chilean Input-Output tables (2021, 12-sector aggregation). The input-output matrix Γi,j

and sectoral material shares (αmi, αvi) come directly from these tables (Subsection 3.5.2).
Consumption parameters (ΩG, ΩS, ωg

i ) are constructed from National Accounts data on

i , ωs

final consumption expenditure by product (Subsection 3.5.3). Home bias parameters (ϱi)

are calculated from trade statistics showing the domestic share of consumption by sector

(Subsection 3.6). Foreign economy steady-state values are normalized: zero net inflation

(Π∗ = 1.00) and the world interest rate is pinned down by the Euler equation (Rw = Π∗/β ≈

1.014 per quarter).

Estimated parameters (Table 3) Several structural parameters are estimated by Sim-

ulated Method of Moments (SMM), matching Chilean aggregate and sectoral business cycle

moments; see Section 3.7 and Appendix F for details. The production input elasticity

19

(ˆεY = 1.484) is close to Cobb-Douglas. The materials-across-sectors elasticity (ˆεm = 0.051)

indicates strong complementarity across IO-linked intermediate input varieties. The import

price adjustment cost (ˆκv = 1,382) approximates near-complete exchange rate pass-through

to import prices. The aggregate labor reallocation cost (ˆc = 0.025) is small, and the ex-

port demand elasticity (ˆη∗ = 0.500) implies inelastic foreign demand for Chilean products.

The goods-services preference shock has moderate persistence (ˆρω = 0.591) and substantial

amplitude (ˆσω = 0.138). TFP shocks share a common persistence of ˆρA = 0.446, with a
cross-sector mean standard deviation of ˆ¯σA = 0.054. The world import price shock has

persistence ˆρP ∗ = 0.515 and standard deviation ˆσP ∗ = 0.056. The foreign demand shock has
ˆρξ = 0.462 and ˆσξ = 0.042. The debt elasticity (ξb = 0.001) and steady-state debt (¯b) are

calibrated rather than estimated.

3.1 Sector-Specific Parameters

Table 4 presents the key production parameters that vary across the 12 sectors of the Chilean

economy.

20

Table 4. Sector-Specific Production Parameters

Sec. Name

αmi

αvi

κi

ϱi

1

2

3

4

5

6

7

8

9

10

11

12

Agriculture, Forestry & Fishing

0.713

0.113

10.47

0.937

Mining

0.745

0.102

4.50

0.999

Manufacturing

0.597

0.282

8.14

0.513

Electricity, Gas, Water & Waste

0.519

0.334

3.60

1.000

Construction

0.608

0.109

1.48

1.000

Trade, Hotels & Restaurants

0.569

0.084

5.14

0.997

Transport, Communications & Info.

0.554

0.163

2.23

0.945

Financial Intermediation

0.568

0.097

1.83

1.000

Real Estate & Housing Services

0.867

0.015

1.25

1.000

Business Services

Personal Services

0.429

0.062

1.97

0.997

0.280

0.048

0.86

1.000

Public Administration

0.228

0.059

0.68

1.000

Notes: αmi is the share of domestic intermediate inputs (materials) in production, αvi is

the share of imported inputs, κi is the Rotemberg price adjustment cost parameter (higher

values indicate greater price stickiness), and ϱi is the home bias parameter in consumption

(Armington weight). Labor share is given by 1 − αmi − αvi. Values are sector-specific and

obtained from data as described in the text. A detailed breakdown of production input

shares is provided in Appendix Table 6.

3.2 Sectoral Consumption Weights

Table 5 shows the consumption weights for goods and services sectors, along with total

consumption levels and shares.

21

Table 5. Sectoral Consumption Weights and Shares

Sec. Name

ωG
i

ωS
i

Ci (Level) Share (%)

1

2

3

4

5

6

7

8

9

10

11

12

Agriculture, Forestry & Fishing

0.073

0.000

0.038

Mining

0.000

0.000

0.000

4.1

0.0

Manufacturing

0.871

0.000

0.392

42.2

Electricity, Gas, Water & Waste

0.056

0.000

0.040

Construction

0.000

0.000

0.000

Retail, Hotels & Restaurants

0.000

0.120

0.056

Transport, Communications & Info.

0.000

0.199

0.082

Financial Intermediation

0.000

0.121

0.058

Real Estate & Housing Services

0.000

0.234

0.104

Business Services

Personal Services

0.000

0.046

0.022

0.000

0.277

0.134

Public Administration

0.000

0.003

0.001

4.3

0.0

6.1

8.9

6.3

11.3

2.4

14.4

0.1

Total

1.000

1.000

0.927

100.0

Notes: ωG

i and ωS

i are consumption weights within goods and services baskets, respec-

tively. Ci represents total consumption of sector i as a fraction of aggregate consumption C.

Share column shows the percentage of total consumption allocated to each sector.

Sectoral consumption weights satisfy P12

i=1 ωG
sified as goods-producing, while sectors 6-12 are service-producing. Weights are computed

i = 1. Sectors 1-5 are clas-

i = 1 and P12

i=1 ωS

from the Chilean National Accounts.

3.3 Production Input Shares

Table 6 provides a detailed breakdown of production input shares across sectors, showing

how each sector combines domestic materials, imported inputs, and labor in its production

function.

22

Table 6. Production Input Shares by Sector - Chile 2021

Sec. Sector Name

αmi

αvi

Labor Total

(Materials)

(Imports) Share Check

1

2

3

4

5

6

7

8

9

10

11

12

Agriculture, Forestry & Fishing

Mining

Manufacturing

Electricity, Gas, Water & Waste

Construction

Trade, Hotels & Restaurants

0.713

0.745

0.597

0.519

0.608

0.569

0.113

0.174

1.000

0.102

0.153

1.000

0.282

0.121

1.000

0.334

0.147

1.000

0.109

0.283

1.000

0.084

0.347

1.000

Transport, Communications & Info.

0.554

0.163

0.283

1.000

Financial Intermediation

Real Estate & Housing Services

Business Services

Personal Services

Public Administration

Mean across sectors

Std. deviation

0.568

0.867

0.429

0.280

0.228

0.557

0.179

0.097

0.335

1.000

0.015

0.118

1.000

0.062

0.509

1.000

0.048

0.672

1.000

0.059

0.713

1.000

0.122

0.321

1.000

0.097

0.193

–

Notes: This table shows the production function input shares for each sector based on the

CES technology: Yi = exp(Ai)

h
mi M (ϵY −1)/ϵY
α1/ϵY

i

+ α1/ϵY

vi V (ϵY −1)/ϵY

i

+ (1 − αmi − αvi)1/ϵY L(ϵY −1)/ϵY

i

iϵY /(ϵY −1)

,

where Mi represents domestic intermediate inputs, Vi represents imported inputs, and Li rep-

resents labor. The labor share is computed as the residual: 1 − αmi − αvi. All values sum

to 1 by construction (shown in ”Total Check” column). Notable patterns: (1) Services sec-

tors (10-12) are relatively labor-intensive with lower material shares; (2) Manufacturing (3)

and Utilities (4) show high import intensity; (3) Real Estate (9) has the highest material

intensity but lowest import dependence. Source: Chilean Input-Output Tables 2021 (MIP

12x12), Central Bank of Chile.

23

Table 7. Production Input Shares by Sector – Chile 2021

Sec. Sector Name

αmi

αvi

αo
vi

α¬o
vi

Labor Total

Materials

Imports

Share Check

1

2

3

4

5

6

7

8

9

10

11

12

Agriculture, Forestry & Fishing

Mining

Manufacturing

Electricity, Gas, Water & Waste

Construction

Trade, Hotels & Restaurants

Transport, Communications & Info.

Financial Intermediation

Real Estate & Housing Services

Business Services

Personal Services

Public Administration

Mean across sectors

Std. deviation

0.713

0.745

0.597

0.519

0.608

0.569

0.554

0.568

0.867

0.429

0.280

0.228

–

–

0.113

0.018

0.095

0.174

0.102

0.022

0.080

0.153

0.282

0.053

0.229

0.121

0.334

0.029

0.305

0.147

0.109

0.005

0.104

0.283

0.084

0.007

0.077

0.347

0.163

0.061

0.102

0.283

0.097

0.000

0.097

0.335

0.015

0.001

0.014

0.118

0.062

0.004

0.058

0.509

0.048

0.002

0.046

0.672

0.059

0.003

0.056

0.713

1

1

1

1

1

1

1

1

1

1

1

1

–

–

–

–

–

–

0.321

1.000

0.193

–

Notes: This table shows the production function input shares for each sector based on the CES

h

technology: Yi = exp(Ai)
,
i
where Mi represents domestic intermediate inputs, Vi represents imported inputs, and Li represents labor.
The labor share is computed as the residual: 1 − αmi − αvi. For each sector i, αvi = αo
vi denoting the

+ (1 − αmi − αvi)1/ϵY L(ϵY −1)/ϵY

mi M (ϵY −1)/ϵY
α1/ϵY

vi V (ϵY −1)/ϵY

vi + α¬o

+ α1/ϵY

i

i

shares of oil and non-oil imports, respectively. All values sum to 1 by construction (shown in ”Total
Check” column). Notable patterns: (1) Services sectors (10-12) are relatively labor-intensive with lower
material shares; (2) Manufacturing (3) and Utilities (4) show high import intensity; (3) Real Estate (9) has
the highest material intensity but lowest import dependence. Source: Chilean Input-Output Tables 2021
(MIP 12x12 and MIP 111x181), Central Bank of Chile.

iϵY /(ϵY −1)

3.4

Input-Output Matrix

Notation convention. We use two complementary representations of the IO structure:

• Model parameter Γi,j: the share of total intermediate inputs of sector i (purchaser)

that are sourced from sector j (supplier). By definition PN

j=1 Γi,j = 1 for all i (each

row sums to 1). This is the convention used in the model equations.

24

• IO table Bi,j: the share of sector i’s (supplier) output absorbed as intermediate inputs

by sector j (purchaser). By definition PN

i=1 Bi,j = 1 for all j (each column sums to 1).

This is the standard input-output table convention, displayed in Table 8.

The two representations are transposes of each other: Γi,j = Bj,i. Table 8 reports Bi,j; the

calibration code transposes this matrix to obtain Γi,j for use in the model.

25

Table 8. Input-Output Matrix (Bi,j, Supplier × Purchaser) – Chile 2021

Producing Sector ↓ / Using Sector →

1

2

3

4

5

6

7

8

9

10

11

12

1. Agriculture

2. Mining

0.001

0.003

0.001

0.078

0.000

0.250

0.305

0.105

0.016

0.217

0.008

0.015

0.002

0.009

0.005

0.005

0.006

0.003

0.003

0.565

0.111

0.008

0.130

0.151

3. Manufacturing

0.006

0.001

0.007

0.002

0.009

0.003

0.003

0.535

0.274

0.002

0.064

0.095

4. Utilities

5. Construction

0.013

0.022

0.202

0.001

0.001

0.099

0.173

0.117

0.016

0.340

0.006

0.009

0.006

0.075

0.003

0.115

0.002

0.001

0.417

0.284

0.048

0.001

0.036

0.012

6. Retail & Hotels

0.187

0.035

0.003

0.529

0.049

0.003

0.004

0.001

0.002

0.003

0.106

0.079

7. Transport & Comm.

0.002

0.014

0.001

0.215

0.015

0.001

0.004

0.309

0.270

0.001

0.108

0.060

2
6

8. Finance

9. Real Estate

10. Business Serv.

11. Personal Serv.

12. Public Admin.

0.001

0.015

0.080

0.058

0.009

0.168

0.411

0.003

0.214

0.001

0.034

0.007

0.000

0.004

0.030

0.108

0.002

0.030

0.039

0.360

0.176

0.242

0.005

0.003

0.002

0.014

0.187

0.057

0.006

0.291

0.000

0.240

0.170

0.001

0.025

0.007

0.023

0.011

0.001

0.160

0.010

0.001

0.306

0.137

0.335

0.001

0.001

0.015

0.014

0.003

0.111

0.169

0.002

0.094

0.167

0.073

0.050

0.287

0.010

0.018

Notes: Element Bi,j represents the share of supplier i’s output used as intermedi-

ate inputs by purchaser j. Each column sums to unity: P12

i=1 Bi,j = 1 for all j. The

model parameter Γi,j = Bj,i (transpose) gives the share of purchaser i’s intermediate inputs
sourced from supplier j; it satisfies P12

j=1 Γi,j = 1 for all i and enters the sectoral price index

(cid:16)P

P M

it =

j Γi,j (P H

jt )1−ϵm
Central Bank of Chile.

(cid:17)1/(1−ϵm)

. Source: Chilean Input-Output Tables 2021 (MIP 12×12),

3.5 Parameters Constructed Directly from Data

Several key structural parameters are obtained directly from Chilean economic data. The

data construction process uses multiple official sources and is documented in the Process

Data Codes folder.

3.5.1 Data Sources and Construction

The model uses data from multiple official sources, processed through scripts in the Process

Data Codes folder. Sectoral GDP, consumption shares, and intermediate input shares from

the National Accounts are processed by a dedicated Stata script. The Input-Output matrix

is constructed and normalized from the MIP 12×12 tables by a Julia preprocessing script

(build IO.jl) and stored in IO 2021 chile.csv. All parameters are loaded at runtime by

the main Julia driver main SOE gap.jl.

The 12 sectors in our model are:

1. Agropecuario-silv´ıcola y Pesca (Agriculture, Forestry and Fishing)

2. Miner´ıa (Mining)

3. Industria manufacturera (Manufacturing)

4. Electricidad, gas, agua y gesti´on de desechos (Utilities)

5. Construcci´on (Construction)

6. Comercio, hoteles y restaurantes (Trade, Hotels and Restaurants)

27

7. Transporte, comunicaciones y servicios de informaci´on (Transport and Communica-

tions)

8. Intermediaci´on financiera (Financial Intermediation)

9. Servicios inmobiliarios y de vivienda (Real Estate)

10. Servicios empresariales (Business Services)

11. Servicios personales (Personal Services)

12. Administraci´on p´ublica (Public Administration)

Data come from the Central Bank of Chile’s Annual National Accounts (Anuario CCNN

2023).1

3.5.2 Input-Output Matrix and Material Shares (Bi,j, Γi,j, αi, αV

i )

The Central Bank of Chile produces comprehensive input-output tables (Matriz Insumo-

Producto, MIP) capturing peso flows between producers and purchasers. We use the 2021

MIP 12×12 table.

The raw MIP table is stored in supplier-by-purchaser format Bi,j (rows = supplying sec-

tors, columns = purchasing sectors; each column sums to 1). To obtain the model parameter

Γi,j = Bj,i (rows = purchasing sectors, each row sums to 1), the matrix is transposed and

row-normalized. The calibration script (IO 2021 chile.csv loaded by main SOE gap.jl)

stores Γi,j directly, so no further transformation is needed at run time.

Concretely, the construction follows three steps:

1. Read the raw 12×12 table Bi,j from mip 12x12.csv (rows = suppliers, columns =

purchasers; each column sums to 1)

2. Transpose: ˜Γi,j = Bj,i (rows = purchasers, columns = suppliers)

3. Normalize each row by its sum to obtain Γi,j, which satisfies P12

j=1 Γi,j = 1 for all i

1https://si3.bcentral.cl/Siete/ES/Siete/Cuadro/CAP_CCNN/MN_CCNN76/CCNN2018_P1/

637801087677220267

28

The resulting Γi,j gives the share of sector i’s intermediate input bundle sourced from sector

j, consistent with the price index P M

it =

(cid:16)P

j Γi,j (P H

jt )1−ϵm

(cid:17)1/(1−ϵm)

in the model.

The material share parameter αi (the share of intermediate inputs in gross production)

is calculated in gen data.do as:

αi =

Intermediatesi
VBPi

where VBP (Valor Bruto de Producci´on) is gross output at producer prices. The labor share

in value added is then 1 − αi.

The import share αV

i represents the share of imported intermediates in total intermediate

consumption by sector i, also obtained from the Input-Output tables.

3.5.3 Consumption Shares

To construct the vectors γg and γs, we use two main data sources: the Production Matrix at

Basic Prices (Table 1) and the Final Use Table at User Prices (Table 6).

We begin by extracting the sectoral shares of goods and services from the production

matrix. These shares are then multiplied by household consumption data from the final use

table to estimate sector-specific expenditures on goods and services.

In this framework:

• Goods correspond to products classified from 1 to 5.

• Services correspond to products classified from 6 to 12.

In the case of Chile, the production matrix is nearly diagonal, indicating that goods are

primarily produced by goods-producing sectors, and services by service-oriented sectors.

The consumption shares within each category are computed as:

γg
i =

Ci
j∈goods Cj

P

,

γs
i =

Ci
j∈services Cj

P

The steady-state weight of goods in the consumption basket, ¯ω, is calibrated to 0.57

based on the relative expenditure shares in the data.

29

3.6 Home Bias in Consumption

The home bias parameter for each sector ϱi is calculated from Chilean national accounts

data contained in the file datos CCNN mayo2025.xlsx. Specifically, ϱi represents the share

of domestically-produced goods in total sectoral consumption (domestic plus imported vari-

eties), computed as:

ϱi =

CH,i
CH,i + CF,i

=

Domestic Consumptioni
Domestic Consumptioni + Imports for Final Consumptioni

where CH,i denotes consumption of the domestic variety and CF,i denotes consumption of

the imported variety of good i.

3.7 Estimation Results

We estimate θ by minimizing the SMM criterion function; see Appendix F for a full de-

scription of the estimator, moment conditions, and identification arguments. Table 9 reports

the estimated values of structural parameters and shock processes. Table 10 compares key
targeted moments against their model counterparts at ˆθ.

Structural parameters. The estimated production input elasticity ˆεY = 1.48 is close to

Cobb-Douglas across labor, domestic materials, and imported inputs. The materials-across-

sectors elasticity ˆεm = 0.051 is well below unity, confirming strong complementarity across

intermediate input varieties — consistent with sectors being largely unable to substitute

away from their IO-determined upstream suppliers in response to relative price changes.

The import price adjustment cost ˆκV ≈ 1,382 effectively approximates the flexible import

prices benchmark used in the impulse response analysis. The aggregate labor reallocation

cost ˆc = 0.025 is small, suggesting that the model requires high measured labor mobility

across sectors to match the observed employment patterns. The export demand elasticity

ˆη∗ = 0.50 implies inelastic foreign demand for Chilean products.

Shock processes. The goods-services reallocation shock is moderately persistent (ˆρω =

0.59) with substantial amplitude (ˆσω = 0.14), indicating that demand shifts between the two

30

consumption aggregates are an important driver of sectoral fluctuations in Chile. Sectoral

TFP shocks have a common persistence of ˆρA = 0.45 and sector-specific standard devia-

tions ranging from 0.002 (Transport, Personal Services) to 0.100 (Construction, Mining),

reflecting significant heterogeneity in supply-side dynamics. The world import price shock

has persistence ˆρP ∗ = 0.52 and standard deviation ˆσP ∗ = 0.056, broadly consistent with the

persistence observed for global commodity prices. The foreign demand shock has ˆρξ = 0.46

and ˆσξ = 0.042.

Moment fit. Table 10 summarizes the fit. The model captures the aggregate volatility of

GDP (0.064 vs. 0.041) and the order of magnitude of inflation volatility (0.039 vs. 0.022). The

aggregate GDP-inflation correlation turns negative in the model (−0.40), consistent with the

sign in the data (−0.02), though substantially more negative, reflecting the dominant role

of supply-side shocks in the estimated model. The model broadly matches services-sector

output volatility but overpredicts goods-sector employment variance. The rank correlations

between model and data sectoral volatilities are near zero across all three margins — out-

put, prices, and employment — indicating that the estimated model does not yet replicate

the cross-sectoral ordering of fluctuations. Improving fit on this dimension is the primary

objective for the next stage of estimation.

4 Results

This section uses the estimated model to analyze the transmission of structural shocks

through the Chilean production network. We study four exercises that illustrate distinct

channels of shock propagation: an aggregate preference shock that shifts the relative de-

mand for goods versus services; a contractionary monetary policy shock; a world oil price

shock that transmits through both the direct import cost channel and the IO network; and

a positive TFP shock to Manufacturing—the most IO-connected goods sector in the Chilean

economy.

For each exercise we compute impulse responses as percentage deviations from the de-

terministic steady state, using the decision rules obtained from the first-order perturbation

31

Table 9. SMM Estimates: Structural Parameters and Shock Processes

Parameter Description

Structural parameters

Estimate

Production input elasticity (labor, materials, imports)

1.484

εY

εm

c

κV

η∗

Materials elasticity (across IO sectors)

Aggregate labor reallocation cost

Import price adjustment cost

Export demand elasticity

Shock processes

ρω

σω

ρA

¯σA

ρP ∗

σP ∗

ρξ

σξ

Preference shock persistence

Preference shock std. dev.

TFP shock persistence (common)

TFP shock std. dev. (cross-sector mean)

World import price shock persistence

World import price shock std. dev.

Foreign demand shock persistence

Foreign demand shock std. dev.

0.051

0.025

1,382

0.500

0.591

0.138

0.446

0.054

0.515

0.056

0.462

0.042

Notes: Parameters estimated jointly by SMM minimizing the weighted distance
between 15 unconditional second moments in Chilean quarterly data (2003Q1–
2023Q4) and their model counterparts; see Appendix F. Bootstrap standard
errors are not yet computed. c = 1/ˆι where ˆι = 40.31 is directly estimated.
κV = exp(ˆℓ) with ˆℓ = 7.23 estimated in logs to enforce positivity.

32

Table 10. SMM Moment Fit: All 45 Targeted Moments

# Sector

Data

Model Data Model Data

Model

Output, std(Yi)

Price, std(P H

i ) Employment, std(Li)

Panel A: Sectoral moments

Goods sectors (1–5)

1 Agriculture & Fishing

0.403

0.057

0.128

0.036

0.125

2 Mining

3 Manufacturing

4 Utilities

5 Construction

Services sectors (6–12)

0.066

0.038

0.063

0.075

0.061

0.167

0.057

0.045

0.051

0.053

0.041

0.030

0.058

0.127

0.059

0.024

0.066

0.035

0.049

0.046

6 Trade & Restaurants

0.075

0.058

0.026

0.036

0.032

7 Transport & Comm.

0.054

0.059

0.031

0.053

0.028

8 Financial Services

9 Real Estate

10 Business Services

11 Personal Services

0.029

0.020

0.060

0.162

0.055

0.047

0.040

0.046

0.069

0.019

0.051

0.029

0.062

0.017

0.064

0.029

0.051

0.101

0.037

0.022

12 Public Administration

0.012

0.052

0.029

0.039

0.047

0.071

0.069

0.093

0.063

0.041

0.040

0.052

0.039

0.041

0.057

0.042

0.043

Panel B: Aggregate moments

std(GDP )

std(π)

corr(GDP, π)

std(Q)

autocorr(Q)

std(T B/GDP )

Data

0.041

0.022

−0.020

0.044

0.716

0.025

Panel C: Cross-sectional rank correlations (Spearman ϱ)

Target

ϱ(cid:0)σY
ϱ(cid:0)σP H
ϱ(cid:0)σL

data, σY
data, σP H
data, σL

model

model

model

(cid:1)

(cid:1)

(cid:1)

1.000

1.000

1.000

Model

0.064

0.039

−0.404

0.035

−0.036

0.059

Model

0.021

−0.077

0.014

Notes: All series HP-filtered (λ = 1600), 2003Q1–2023Q4. Panel A reports sector-by-sector standard
deviations. Panel B reports scalar aggregate moments. Panel C reports Spearman rank correlations
between the 12-sector vectors of model and data standard deviations; a target of 1 indicates perfect
cross-sectoral ordering. Estimated parameters ˆθ reported in Table 9; SMM objective value at ˆθ:
55.2.

33

of the model. Shock sizes are normalized to facilitate comparison with empirical estimates:

one standard deviation for the preference shock (σξ), a 25 basis point unexpected tightening

for the monetary policy shock, a 10% world oil price increase, and a 1% TFP improvement

for the Manufacturing shock. All impulse responses are computed using the Julia script

oil shock analysis.jl.

4.1 Aggregate Preference Shock

A positive preference shock (ξt > 1) raises the marginal utility of current consumption, ex-

panding aggregate demand uniformly across sectors. Because price stickiness varies across

sectors, the real allocation effects differ: sectors with stickier prices (higher κi, such as Agri-

culture and Manufacturing) expand output more on impact, while services sectors—which

adjust prices more quickly—accommodate a larger share of the demand increase through

prices rather than quantities. The IO network transmits the demand expansion upstream:

as final goods sectors expand, they demand more intermediate inputs from their suppliers,

generating a second round of output and price increases that propagates throughout the

network. The monetary authority responds by tightening the policy rate, which partially

crowds out the demand expansion and contributes to a real exchange rate appreciation that

dampens export demand.

4.2 Monetary Policy Shock

A contractionary monetary policy shock (an unexpected increase in the policy rate Rt)

tightens financial conditions and reduces consumption through the intertemporal substitu-

tion channel. The decline in aggregate demand reduces output across all sectors, but the

adjustment is uneven: sectors with stickier prices experience larger output contractions (as

they cannot lower prices quickly to restore demand), while sectors with more flexible prices

see a larger price decline. The real exchange rate appreciates as the domestic interest rate

rises relative to the world rate, consistent with the UIP condition. The appreciation lowers

the domestic-currency cost of imports, which partially offsets the contractionary demand

effect on inflation. The trade balance deteriorates as export demand falls in response to the

34

exchange rate appreciation.

4.3 Oil Price Shock

We study the transmission of a world oil price shock through the Chilean economy. The

exercise is motivated by Chile’s dependence on imported energy: crude oil and refined fuels

constitute a significant share of intermediate imports across most sectors, with particu-

larly high exposure in Transport and Communications (37%), Mining (22%), Manufacturing

(19%), and Agriculture (16%). Because oil enters production as an intermediate input, the

first-round cost increase in oil-intensive sectors propagates to their customers through the

input-output (IO) network—a general-equilibrium amplification channel that is absent in

representative-agent or single-sector models. We quantify this channel by decomposing the

full model response into a direct cost effect and a network amplification effect, and connect

both components to closed-form approximations that we derive from the model’s marginal

cost equations.

4.3.1 Experiment Design

We subject the economy to a one-time 10% increase in the world oil price P O∗

t

, governed by

the AR(1) process

log

!

  P O∗
t
P O∗

= ρO∗ log

!

  P O∗
t−1
P O∗

+ σO∗ εO∗
t

,

(26)

with persistence ρO∗ = 0.9 and innovation standard deviation σO∗ = 0.02. All other shocks

are set to zero. The domestic oil price satisfies P O

t = EtP O∗

t

, so that movements in the

nominal exchange rate fully transmit the foreign-currency price change to domestic costs.

In each sector i, composite intermediate imports Vit are a CES aggregate of oil (V Oil

it ) and

non-oil (V Non

it

) imports with elasticity of substitution ϵoil

V = 0.5, reflecting limited short-run

substitutability between energy and other imported inputs. The associated import price

index satisfies

(cid:16)

P V
it

(cid:17)1−ϵoil

V = αOil

i

(cid:16)

P O
t

(cid:17)1−ϵoil

V + (1 − αOil
i )

(cid:16)

P v
t

(cid:17)1−ϵoil
V ,

(27)

35

where the sector-specific oil shares αOil

i

are drawn from Chile’s 2021 Input-Output tables.2

4.3.2 Transmission Channels and Derivation of the Approxima-

tions

The oil price shock transmits to the economy through two distinct channels. Before charac-

terizing the quantitative results, we derive closed-form expressions for each channel directly

from the model’s production block. These expressions clarify the sources of heterogeneity in

sectoral responses and connect the model to the Leontief inverse framework of ?.

Step 1: Log-linearizing marginal cost. The marginal cost of a representative interme-

diate firm in sector i is

M C i

t =

(cid:18)
αmi (P M i

t

1
Ai
t

)1−ϵY + αV i (P V

it )1−ϵY + (1 − αmi − αV i) (P Li

t )1−ϵY

(cid:19) 1
1−ϵY ,

(28)

where P M i

t

, P V

it , and P Li

t are the prices of the materials bundle, composite imported inputs,

and labor services, respectively. Log-linearizing (28) around the deterministic steady state

— where all prices are normalized to unity so that the CES cost-share weights coincide with

the calibrated factor shares — yields

dmci

t = αmi bpM i

t + αV i bpV

it + (1 − αmi − αV i) bpLi

t − ˆAi
t,

(29)

where bxt ≡ log xt − log ¯x denotes the log deviation from steady state.

Step 2: Log-linearizing the oil import price index. Log-linearizing (27) around the

steady state gives

it = αOil
bpV

i

t + (1 − αOil
bpO

i ) bpv
t .

(30)

When the only shock is the world oil price ( bpv
i is simply αOil

i

times the oil price increase, scaled by the imported-input share αV i in the

t = 0), the import price response in sector

marginal cost equation.

2Specifically, we use the use matrix (Cuadro 21) for products 28 and 77–81 (crude oil, diesel, gasoline,

kerosene, fuel oils, LPG) at basic prices, aggregated to our 12-sector classification.

36

Step 3: Direct cost channel. To isolate the first-round effect, suppose the prices of all

t ), non-oil imports ( bpv
other inputs — labor ( bpLi
yet adjusted, and abstract from TFP shocks ( ˆAi

t ), and the materials bundle ( bpM i

t

) — have not

t = 0). Substituting (30) into (29) then gives

∆ log M C i
t

(cid:12)
(cid:12)
(cid:12)
(cid:12)direct

≈ αV i αOil

i ∆ log P O
t .

(31)

The direct cost increase is therefore the product of two sector-specific parameters: the share

of production cost attributable to imported inputs (αV i) and the oil intensity of those imports

(αOil

i ). This product varies by a factor of roughly sixty across our twelve sectors — from 0.003

for Financial Services to 0.192 for Transport and Communications — generating substantial

heterogeneity in the first-round inflationary pressure.

Step 4: Network amplification channel. The materials bundle of sector i is a CES

aggregate of the domestic output of all N sectors:

P M i

t =





N
X

j=1



1
1−ϵm

Γi,j (P H

jt )1−ϵm



,

(32)

where Γi,j is the share of sector i’s material inputs purchased from sector j, with P
Log-linearizing at the steady state (where P H

j = 1 for all j) yields

j Γi,j = 1.

bpM i
t =

N
X

j=1

Γi,j bpH

jt =

(cid:16)

Γ bpH

t

(cid:17)

,

i

(33)

so the materials price of sector i is a weighted average of the output prices of all sectors that

supply it.

To close the system, we use the observation that in the steady state the Rotemberg pricing

problem implies a constant markup over marginal cost: P H

relationship gives bpH
point of approximation.3 Substituting bpH

ϵ−1 M C i. Log-linearizing this
i = dmci, so output prices move one-for-one with marginal costs at the
t into (33) and then into (29), and isolating

jt = dmcj

i = ϵ

3This static pass-through condition holds exactly at the zero-inflation deterministic steady state. The full
dynamic model incorporates the inertial Rotemberg adjustment costs κi, which cause partial and delayed
pass-through in general; the approximation therefore captures the impact-period response in the limiting
case of flexible prices.

37

the oil shock ( bpv

t = bpLi

t = ˆAi

t = 0), we obtain the recursive marginal cost system:

dmci

t = αmi

N
X

j=1

Γi,j dmcj

t + αV i αOil

i

bpO
t ,

i = 1, . . . , N.

(34)

Equation (34) makes the network propagation transparent: the marginal cost of sector i

rises both because oil directly raises its own import costs (the second term) and because the

cost increases in its upstream suppliers are passed forward through the IO matrix (the first

(cid:12)
(cid:12)
term). Writing (34) in vector form, dmct = diag(αm) Γ dmct + dmc
(cid:12)direct

, and solving gives

(cid:12)
(cid:12)
∆ log MC
(cid:12)
(cid:12)total

(cid:16)

=

I − diag(αm) Γ

(cid:17)−1

(cid:12)
(cid:12)
∆ log MC
(cid:12)
(cid:12)direct

,

(35)

where the Leontief inverse (I − diag(αm) Γ)−1 encodes all rounds of network propagation.

The (i, j) element of this inverse measures the total increase in sector i’s marginal cost per

unit of direct cost shock to sector j, accounting for all indirect IO linkages. The network

amplification effect for sector i is the difference between the total and direct effects:

∆ log M C i

(cid:12)
(cid:12)
(cid:12)
(cid:12)network

= ∆ log M C i

(cid:12)
(cid:12)
(cid:12)
(cid:12)total

− ∆ log M C i

(cid:12)
(cid:12)
(cid:12)
(cid:12)direct

.

(36)

Two conditions guarantee that the Leontief inverse (I −diag(αm) Γ)−1 exists and has non-

negative entries. First, since αmi < 1 for all i (intermediate materials do not exhaust the

entire production budget) and Γ is a row-stochastic matrix, the spectral radius of diag(αm) Γ

is strictly less than one, ensuring invertibility. Second, both diag(αm) and Γ are non-negative,

so by the Perron-Frobenius theorem the inverse is also non-negative: cost shocks in one sector

cannot reduce costs in another. We note that the approximation relies on three maintained

assumptions: (i) prices equal markups over marginal cost (static pass-through); (ii) the

non-oil import price and wages are fixed; and (iii) TFP is unaffected by the shock. The

full dynamic model relaxes all three: Rotemberg adjustment costs introduce sluggish pass-

through, the exchange rate endogenously affects non-oil import prices, and the labor market

equilibrium is determined jointly with production. The approximation therefore serves as a

transparent diagnostic tool rather than an exact decomposition of the model IRFs.

38

4.3.3 Aggregate Responses

Figure 1 reports the impulse responses of aggregate variables to the 10% oil price shock. GDP

falls by 0.19% on impact and remains below trend for several quarters, recovering gradually

as the shock dissipates. CPI inflation rises by 0.15 annualized percentage points on impact,

peaking at 0.21 pp in the second quarter, as higher import costs pass through to domestic

prices via both the direct channel (31) and the IO multiplier (35). The central bank responds

by raising the policy rate by 0.10 annualized percentage points on impact, peaking at 0.23

pp by the third quarter; the resulting real appreciation — Q falls by 0.50% on impact

— partially offsets the domestic-currency oil price increase. The trade balance improves

markedly (0.44% of GDP on impact, peaking at 2.40% of GDP by the third quarter), as

demand compression reduces import volumes by more than the higher oil price raises the

import bill.

Table 11 summarizes the impact and peak responses of the key aggregate variables.

All impulse response figures in this section overlay three parameterizations of the production-

input substitution elasticity ϵY : a near-Leontief specification (ϵY = 0.30), the baseline cal-

ibration (ϵY = 0.80), and a high-substitution case (ϵY = 2.0). This comparative statics

exercise isolates a key structural channel: when ϵY is low, firms cannot substitute away

from costlier oil-intensive input bundles, so the marginal cost increase is large and persists

across sectors; when ϵY is high, firms reallocate toward cheaper inputs, compressing the cost

pass-through and attenuating the aggregate stagflationary response. As Figure 1 shows, the

near-Leontief case generates a substantially deeper GDP contraction and a larger inflation-

ary spike, while the high-substitution case moderates both. The labor market responses and

the cross-sectoral dispersion in output and inflation are similarly sensitive to this parameter,

as we document below.

Figure 2 aggregates the sectoral home-price inflation responses into goods, services, and

economy-wide averages, weighted by steady-state consumption expenditure shares. Goods

home-price inflation reaches 1.19 annualized percentage points on impact, driven by the direct

oil cost channel operating through Transport, Manufacturing, and Mining. Services home-

price inflation is more muted at 0.33 pp on impact, as services sectors have lower direct oil

39

Figure 1. Aggregate impulse responses to a 10% world oil price shock. Each panel shows three
lines corresponding to different values of the production-input substitution elasticity: near-
Leontief (ϵY = 0.30, green dashed), baseline (ϵY = 0.80, blue solid), and high substitution
(ϵY = 2.0, red dotted). All variables in percentage deviations from steady state. The shock
follows an AR(1) with persistence ρO∗ = 0.9.

Table 11. Respuestas Agregadas a un Shock de 10% al Precio del Petr´oleo

Variable

Impacto M´ın/M´ax Trimestre

PIB

-0.491

-0.491

Inflaci´on IPC

+0.070

+0.106

Consumo

Empleo

-0.256

-0.264

+0.095

+0.095

Tipo de Cambio Real

-0.103

+0.147

Balanza Comercial

-0.240

-0.240

Tasa de Pol´ıtica

+0.045

+0.139

1

2

2

1

37

1

5

Notas: PIB, consumo, empleo y TCR en % de desviaci´on del
estado estacionario; inflaci´on y tasa de pol´ıtica en pp anualizados;
balanza comercial en pp del PIB. El shock es un aumento ´unico
de 10% en el precio mundial del petr´oleo P O∗
con persistencia
Impacto = respuesta en el trimestre 1. M´ın/M´ax =
ρ = 0.9.
respuesta extrema en 40 trimestres.

t

40

intensity and the IO spillovers from goods to services are dampened by greater price stickiness

in several services sectors. The expenditure-weighted aggregate home-price inflation is 0.82

pp on impact — roughly four times larger than the aggregate CPI response. The wedge

reflects a key open-economy offset: the real appreciation induced by monetary tightening

lowers the consumer price of the foreign variety in every Armington bundle, substantially

insulating consumers from the domestic cost-push shock even as home-price inflation across

sectors remains substantial.

Figure 2. Expenditure-weighted home-price inflation by sector group. All series in annual-
ized pp deviations from steady state. Aggregate (solid black): economy-wide consumption-
share-weighted average. Goods (dashed blue): within-group weighted average, sectors 1–5.
Services (dotted red): within-group weighted average, sectors 6–12.

4.3.4 Sectoral Responses

The heterogeneous oil exposure across sectors generates a rich pattern of cross-sectoral re-

sponses. Figure 3 displays the output IRFs for all 12 sectors. The sectors most affected are

those with the highest direct oil intensity — Transport and Communications, Mining, and

Manufacturing — which face the largest first-round cost increases. However, the ranking of

output contractions does not mirror the ranking of oil shares: sectors with low direct oil ex-

posure but high reliance on intermediate inputs from oil-intensive upstream sectors — such

as Trade and Business Services — experience significant contractions through the network

channel in (35).

41

Figure 3. Sectoral output responses to a 10% oil price shock. Each subplot shows three
lines for near-Leontief (ϵY = 0.30, green dashed), baseline (ϵY = 0.80, blue solid), and high
substitution (ϵY = 2.0, red dotted).

42

Figure 4 shows the sectoral price level responses. Sectors with higher Rotemberg ad-

justment costs κi — stickier prices — adjust more slowly, absorbing the cost shock through

compressed markups in the short run rather than immediate price increases. The hetero-

geneity in κi (see Section 3), combined with the IO structure, implies that the inflationary

impulse propagates unevenly across sectors even when the underlying cost shock is common

in the input space.

Figure 4. Sectoral home-price level responses to a 10% oil price shock. All series in percentage
deviations from steady state.

Figure 5 translates the price-level responses into sectoral nominal home-price inflation

rates. Sectoral inflation is defined as the gross quarter-over-quarter change in the nominal

home price of sector i:

ˆΠH

it ≡ Πt ·

P H
it
P H

it−1

,

(37)

where P H
it

is the relative home price (normalized by aggregate CPI Pt) and Πt = Pt/Pt−1

43

is gross CPI inflation. Reported IRFs are annualized quarter-over-quarter log changes
(cid:16)
4

. The cross-sectoral dispersion is substantial. Oil-intensive goods sec-

∆ log P H

(cid:17)

it + log Πt

tors record the largest inflation responses on impact: Transport and Communications leads,

followed by Manufacturing and Mining, all reflecting high direct oil shares. Agriculture and

Utilities also see meaningful inflation despite lower direct oil exposure, because they source

heavily from oil-intensive upstream suppliers.

In contrast, several services sectors record

negative inflation on impact: Financial Services, Personal Services, Public Administration,

and Business Services all fall below zero. These sectors have minimal direct oil exposure and

are subject to demand compression from the monetary tightening; with relatively lower price

stickiness, they adjust downward quickly through the demand-compression channel. The sign

reversals illustrate the two counteracting forces in oil shock transmission: cost-push domi-

nates for oil-intensive sectors, while demand-compression dominates for oil-remote services

sectors.

4.3.5 Labor Market Responses

The oil shock generates a contraction in aggregate labor demand alongside a decline in real

wages. Figure 6 reports the impulse responses of aggregate employment Nt and the real wage

wt. Employment falls on impact and remains below trend for several quarters, mirroring the

GDP contraction: as higher marginal costs compress production across sectors, firms reduce

their labor input in proportion to the output decline. The real wage also declines on impact,

reflecting the fall in labor’s marginal revenue product as production costs rise and output

contracts. The wage response is somewhat more persistent than the employment response,

as the labor market clears more slowly through wage adjustment than through quantity

adjustment under the calibrated labor supply elasticity.

The sectoral employment responses, displayed in Figure 7, reveal substantial cross-

sectoral heterogeneity that reflects the interaction of the cost-push channel with sectoral

labor demand elasticities and the IO network. Oil-intensive sectors—Transport and Com-

munications, Mining, and Manufacturing—record the largest employment contractions on

impact, consistent with their severe marginal cost increases. However, the sectoral employ-

ment ranking does not perfectly mirror the output ranking: sectors with high labor intensity

44

Figure 5. Sectoral nominal home-price inflation (annualized pp) responses to a 10% world oil
price shock. Inflation is defined as ˆΠH
it−1 (equation (37)), reported as annualized
quarter-over-quarter log changes 4(∆ log P H
it + log Πt). Blue panels: goods sectors (1–5); red
panels: services sectors (6–12).

it = Πt·P H

it /P H

Figure 6. Aggregate labor market responses to a 10% oil price shock. Left panel: aggregate
employment Nt; right panel: real wage wt. Three ϵY cases overlaid as in Figure 1.

45

(low αmi) but low direct oil exposure may experience relatively larger employment contrac-

tions than their output loss would suggest, because the labor input absorbs a disproportionate

share of the cost adjustment when material and import inputs are relatively complementary

(ϵY < 1). Conversely, material-intensive sectors with high oil exposure substitute away from

the costlier material bundle toward labor, partially cushioning sectoral employment even as

output contracts. The sensitivity to ϵY is particularly pronounced in the labor market: under

near-Leontief technology, firms cannot substitute oil for labor, so the employment contrac-

tion in oil-intensive sectors is amplified; under high substitution, firms shift toward labor as

oil becomes costlier, partially cushioning the employment decline even as output contracts.

The cross-sectoral dispersion in employment responses is economically meaningful: the most

affected sector contracts by several times the economy-wide average, underscoring that ag-

gregate labor market statistics can mask substantial sectoral reallocation pressures induced

by commodity price shocks.

4.3.6 Output Gap Responses

The output gap—the deviation of actual output from its flexible-price counterpart—isolates

the component of the output decline that is attributable to nominal rigidities and monetary

policy, net of the efficient response to the real shock. Figure 8 reports the aggregate output

gap and GDP gap. The oil shock opens a negative output gap: sticky prices prevent the

economy from adjusting to the new efficient allocation, and the Taylor rule’s response to the

inflationary impulse further tightens monetary conditions, amplifying the real contraction

beyond what would obtain under flexible prices. The gap is largest under the near-Leontief

specification (ϵY = 0.30), where the inability to substitute away from oil-intensive inputs

magnifies both the inflationary pressure—triggering a stronger monetary tightening—and

the real cost of price stickiness. Under high substitution (ϵY = 2.0), the output gap is

substantially smaller, as firms can reallocate across inputs more easily, compressing both the

cost pass-through and the monetary response.

The sectoral output gaps, shown in Figure 9, reveal that the welfare cost of nominal

rigidity is unevenly distributed across sectors. Sectors with high price stickiness (κi) exhibit

larger output gaps because their prices adjust sluggishly, distorting relative prices and am-

46

Figure 7. Sectoral employment responses to a 10% oil price shock. Three ϵY cases overlaid
as in Figure 1.

47

Figure 8. Aggregate output gap and GDP gap responses to a 10% oil price shock. Three ϵY
cases overlaid as in Figure 1.

plifying misallocation. The interaction between sectoral price stickiness and the IO network

creates an additional channel: a sector whose upstream suppliers have sticky prices faces a

distorted input cost signal, which propagates through the production network and widens

the downstream output gap even when the downstream sector’s own prices are relatively

flexible. This pattern is especially visible in services sectors such as Business Services and

Personal Services, which have moderate direct oil exposure but large output gaps driven by

upstream price distortions transmitted through the IO network.

4.3.7 Decomposition: Direct Cost vs. Network Amplification

Table 12 decomposes the marginal cost increase in each sector into the direct channel (31)

and the network amplification (36), as computed from equations (31)–(35). Two key findings

emerge.

First, the network amplification is quantitatively large. On an output-weighted basis,

the total marginal cost increase is 2.50 times the direct oil cost effect: the weighted-average

direct effect is 0.24% and the total is 0.61%, so IO propagation accounts for approximately

48

Figure 9. Sectoral output gap responses to a 10% oil price shock. Three ϵY cases overlaid as
in Figure 1.

49

60% of the aggregate inflationary impact. Ignoring input-output linkages would therefore

cause a researcher to underestimate the inflationary impulse of an oil shock by a factor of

2.5.

Second, the network channel partially reshuffles the cross-sectoral vulnerability ranking.

While Transport and Communications has the highest direct exposure, Agriculture and

Manufacturing receive substantial additional cost pressure through the IO network because

they purchase heavily from other oil-intensive sectors. The amplification ratio is highest for

sectors with low direct oil exposure but high upstream dependence on oil-intensive suppliers

— precisely the sectors for which the Leontief inverse captures dynamics that the direct

channel alone would miss.

Figure 10 provides a complete picture of the two-stage transmission from the world

oil price to domestic consumer prices. The stacked bars decompose the total marginal

cost increase in each sector into the direct oil cost channel (31) (orange) and the network

amplification channel (36) (blue); their combined height equals the model-based marginal

cost response from the Dynare general-equilibrium solution, expressed as a percentage de-

viation from steady state. The black diamonds report the nominal home-price inflation
log ˆΠH

i,1 + log Π1 (equation (37)) at the impact quarter, in quarter-over-quarter

i,1 = ∆ log P H

percentage points, also from the Dynare IRF.

The first thing the figure reveals is the cross-sectoral pattern of cost exposure. Transport

and Communications faces the largest total marginal cost increase, driven primarily by its

high direct oil share; oil-remote sectors such as Financial Services and Business Services

face negligible or negative total cost effects, the latter reflecting demand compression from

monetary tightening that reduces variable input costs. The network amplification (blue) is

sizable in every oil-intensive sector: Agriculture and Manufacturing, for instance, receive

substantial additional cost pressure through upstream IO linkages on top of their direct

exposure, consistent with the factor-of-2.5 aggregate amplification discussed above.

The second thing the figure reveals is the pass-through gap between the bar tops and the

diamonds. This gap is not a model artifact: it reflects the Rotemberg price adjustment costs

κi, which penalize rapid price increases and therefore cause firms to partially absorb cost

shocks through markup compression rather than immediate pass-through. Log-linearizing

50

the Rotemberg pricing condition around the zero-inflation steady state gives the sectoral

NKPC:

log ˆΠH

i,t ≈

ε − 1
κi

dmci

t + β Et

h

log ˆΠH

i,t+1

i

,

(38)

where ˆΠH

it = Πt · P H

it /P H

it−1 (equation (37)). The pass-through at impact is proportional

to (ε − 1)/κi: a high κi compresses the contemporaneous coefficient, producing a large

bar-to-diamond gap; a low κi allows rapid adjustment, closing the gap. Agriculture (κ1 =

10.47) and Manufacturing (κ3 = 8.14) have among the highest adjustment costs, so their

diamonds fall well short of their bar tops. Transport and Communications (κ7 = 2.23) has

substantially lower adjustment costs; its diamond sits close to the bar top, indicating near-

complete impact-period pass-through. Sticky services sectors such as Financial Services and

Business Services show negative diamonds below bars that already sit near zero: monetary

tightening compresses demand for these sectors even as cost pressure from oil is minimal,

pushing their inflation below the pre-shock level.

Figure 11 plots the impact-period sectoral output and marginal cost responses side by

side. The juxtaposition reveals two patterns that the aggregate numbers mask. First, the

sectors with the largest marginal cost increases — Transport and Communications, Mining,

and Manufacturing — are not necessarily those with the largest output contractions: sectoral

output depends on the balance between the cost-push effect and the endogenous demand

response, which in turn reflects the general-equilibrium adjustment of wages, the exchange

rate, and aggregate consumption. Second, for several services sectors the marginal cost

increases modestly while output contracts noticeably, reflecting their exposure through the

demand channel rather than the cost channel.

4.3.8 The Role of Sectoral Heterogeneity

The quantitative importance of heterogeneity — in oil exposure, price stickiness, and net-

work centrality — motivates the multi-sector framework. A representative-sector model that

uses economy-wide averages of αOil and κ would miss both the cross-sectoral dispersion in

inflationary pressure and the sign reversals documented in Figure 5. Equally, a multi-sector

model without IO linkages would capture the cross-sectoral heterogeneity in first-round costs

51

(a) Near-Leontief (ϵY = 0.30)

(b) Baseline (ϵY = 0.80)

52

(c) High substitution (ϵY = 2.0)

Figure 10. GE marginal cost decomposition and inflation pass-through at impact (t = 1)

following a 10% oil price shock, under three values of the production-input substitution

elasticity ϵY . Orange: direct oil cost channel. Blue: network propagation via intermediate

prices. Green: non-oil imports. Crimson:

labor cost. Gray:

residual (GE effects not

captured by the partial decomposition). Black diamonds: nominal home-price inflation at

impact. Under near-Leontief technology (panel a), firms cannot substitute away from oil,

so the direct cost channel dominates and marginal cost increases are large. Under high

substitution (panel c), firms reallocate toward non-oil inputs, compressing the cost pass-

through and reducing the inflationary impulse.

Table 12. Descomposici´on Sectorial de un Shock de 10% al Precio del Petr´oleo

Sector

Part. (%) Directo Red Total Raz´on

IRF (%)

Petr´oleo

Aumento CM (%)

Amplif. Producto

Agricultura y Pesca (G)

Miner´ıa (G)

Manufactura (G)

Electricidad (G)

Construcci´on (G)

Comercio y Hoteles (S)

16.4

21.6

18.9

8.7

4.2

7.9

0.176

0.504

0.680

3.87

+0.080

0.210

0.472

0.682

3.24

+0.110

0.509

0.469

0.978

1.92

+0.032

0.277

0.320

0.598

2.16

+0.074

0.043

0.456

0.500

11.58

−0.114

0.064

0.368

0.432

6.77

−0.043

Transporte y Com. (S)

37.3

0.582

0.422

1.004

1.72

+0.143

Finanzas (S)

Inmobiliario (S)

Serv. Empresariales (S)

Serv. Personales (S)

Adm. P´ublica (S)

Agregado (ponderado por producto)

0.4

4.5

6.3

3.9

4.5

—

0.004

0.173

0.177

44.67

−0.086

0.006

0.387

0.393

62.13

−0.120

0.037

0.197

0.234

6.30

+0.011

0.018

0.176

0.193

10.89

−0.168

0.025

0.138

0.163

6.47

−0.093

0.286

0.395

0.680

2.38

−0.491

Notas: La participaci´on del petr´oleo es αOil
, la fracci´on de las importaciones intermedias del sector i que
i ∆ log P O. La
corresponde a petr´oleo/combustibles (matriz IP de Chile 2021). El CM directo es αV i αOil
red es el costo adicional v´ıa encadenamientos IP bajo la aproximaci´on de Leontief de equilibrio parcial.
Total = Directo + Red. Raz´on de amplificaci´on = Total/Directo. La IRF de producto es la respuesta
sectorial en el per´ıodo de impacto. B = Bienes, S = Servicios.

i

53

(a) Near-Leontief (ϵY = 0.30)

(b) Baseline (ϵY = 0.80)

(c) High substitution (ϵY = 2.0)

Figure 11. Impact-period sectoral output and marginal cost responses to a 10% oil price
54
shock under three values of ϵY . Blue bars: output Yi; red bars: marginal cost M Ci; both in
percentage deviations from steady state at t = 1. With low substitution (panel a), output
contractions are larger and more dispersed because firms cannot substitute away from costlier
oil-intensive inputs. With high substitution (panel c), the output response is attenuated as

firms reallocate toward cheaper input bundles.

but miss the network amplification that accounts for the majority of the aggregate inflation-

ary impulse.

Figure 12 documents the oil import shares across sectors. The dispersion is large: Trans-

port and Communications allocates 37% of its imported inputs to oil, while Financial Services

allocates less than 1%. This variation, combined with the low substitutability between oil

and non-oil imports (ϵoil

V = 0.5 < 1), implies that oil-intensive sectors cannot easily substitute

away from oil in the short run, making the cost pass-through approximately proportional to

the direct oil share.

Figure 12. Oil share of intermediate imports by sector (αOil
IO tables.

i ), computed from Chile’s 2021

The network amplification depends on the density and topology of the IO matrix. Sectors

that are central in the production network — those supplying intermediate inputs to many

other sectors — amplify the oil shock even if their own direct oil exposure is moderate.

Conversely, peripheral sectors experience effects close to the direct channel alone. This

echoes the theoretical results in ? on the role of network centrality in shock propagation,

and here we quantify the channel in a calibrated model with realistic Chilean IO linkages

and heterogeneous price rigidities.

55

4.4 Manufacturing TFP Shock

We now study the transmission of a positive technology shock originating in the Manu-

facturing sector (sector 3). This experiment serves two purposes. First, it isolates the

demand-pull and cost-relief spillover channels that the IO network generates in response to

a supply shock originating within the domestic economy—as opposed to the external cost-

push channel studied in Section 4.3. Second, Manufacturing is the most IO-connected goods

sector in Chile’s production network, supplying intermediate inputs to all other goods sec-

tors and several services sectors; a productivity improvement here therefore exercises the full

upstream-to-downstream propagation machinery of the model.

4.4.1 Experiment Design

We impose a one-time 1% increase in the TFP level of the Manufacturing sector, governed

by the AR(1) process

log

!

  A3
t
A3

= ρA log

!

  A3
t−1
A3

+ σA,3 εA3
t

,

(39)

with estimated persistence ρA = 0.446 and sectoral standard deviation σA,3 = 0.099, both

from the SMM estimation. All other sectoral TFP processes and external shocks are held

at zero. The 1% shock magnitude is chosen to facilitate direct comparison with the 10% oil

price shock on a normalized basis.

4.4.2 Transmission Channels

A positive TFP shock in Manufacturing transmits to the rest of the economy through three

distinct channels.

Direct supply channel. Higher TFP in sector 3 shifts out the production possibility

frontier of Manufacturing, lowering its real marginal cost:

∆ log M C 3
t

(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)direct

≈ −∆ log A3

t < 0.

(40)

56

Under Rotemberg pricing, the lower marginal cost feeds into lower sectoral home-price in-

flation ΠH

3,t, which acts as a domestic deflationary force. The output of sector 3 expands

both because productivity raises supply capacity and because the lower relative price stim-

ulates demand from households, exporters, and other sectors purchasing Manufacturing as

an intermediate.

Downstream cost-relief channel. The fall in Manufacturing home prices P H

3,t reduces

the intermediate input cost of every sector that purchases Manufacturing as a material

input, according to the IO matrix Γi,3. The first-order approximation of the pass-through to

downstream sector i is

∆ log P M i

t ≈ αmi Γi,3 ∆ log P H
3,t,

(41)

where αmi is the share of materials in sector i’s production costs and Γi,3 is the purchasing

weight on Manufacturing inputs. Sectors with high material intensity and strong Manufac-

turing linkages—Construction, Utilities, and Agriculture—experience the largest cost-relief

spillovers. Including second-round network propagation via the Leontief inverse, the total

downstream deflation is larger than the first-order expression above.

Demand-spillover channel. The expansion in Manufacturing output raises household

income and labor demand in sector 3. Through the aggregate labor market and the household

budget constraint, higher income stimulates consumption demand across all sectors. This

demand effect is partially inflationary for services sectors, where the cost-relief channel is

weaker (services purchase relatively little Manufacturing input) but the demand channel

reaches them through final consumption. The net price response in services therefore depends

on the balance between these two opposing forces.

4.4.3 Aggregate Responses

Figure 13 reports the impulse responses of aggregate variables to the 1% Manufacturing TFP

shock. The macroeconomic transmission is expansionary-deflationary—the mirror image of

the stagflationary oil price shock.

57

GDP rises by 0.39% on impact as higher Manufacturing productivity directly expands

output in the sector-of-origin and stimulates downstream production through the IO network.

The response is hump-shaped, with GDP remaining above trend for several quarters before

mean-reverting under the estimated AR(1) persistence of ρA = 0.446.

Consumer price inflation falls by 0.39 annualized percentage points on impact, driven

by the direct pass-through of lower Manufacturing prices into the consumer basket and by

the downstream cost relief that propagates deflationary pressure to upstream goods sectors.

The central bank responds by easing: the policy rate falls by 0.25 annualized percentage

points on impact, reaching its trough in the second quarter. This monetary accommodation

induces a real exchange rate depreciation—Q rises by 0.19% on impact—which partially

offsets the domestic deflationary impulse through higher import prices, but leaves the net

inflation response firmly negative.

The trade balance improves by 0.28% of GDP on impact: the productivity gain raises

export competitiveness directly (lower unit costs for Manufacturing exports) and indirectly

through the exchange rate depreciation, more than offsetting the modest increase in import

volumes induced by the income expansion.

Figure 13. Aggregate impulse responses to a 1% Manufacturing TFP improvement. All
variables in percentage deviations from steady state (GDP, Q, TB/GDP) or annualized
percentage-point deviations (CPI inflation, policy rate). The shock follows an AR(1) with
persistence ρA = 0.446.

58

4.4.4 Sectoral Responses

The sectoral pattern of output responses, shown in Figure 14, reflects the interplay of direct

supply expansion, downstream cost relief, and aggregate demand spillovers. Manufactur-

ing output itself rises sharply on impact, consistent with the direct productivity increase.

Goods sectors that use Manufacturing as a primary intermediate—Construction, Utilities,

and Agriculture—expand meaningfully through the downstream cost-relief channel:

lower

Manufacturing input prices reduce their marginal costs, stimulating production at given de-

mand. Services sectors also expand on impact, driven primarily by the aggregate demand

spillover as rising household income shifts out consumption demand for all varieties.

Figure 14. Sectoral output responses to a 1% Manufacturing TFP shock. Blue panels: goods
sectors (1–5); red panels: services sectors (6–12).

Figure 15 displays the sectoral nominal home-price inflation responses, defined as ˆΠH

it =

Πt · P H

it /P H

it−1 (equation (37)). The cross-sectoral pattern reveals the contrasting force of

59

supply and demand channels. Manufacturing itself records the largest nominal price decline

on impact, as the direct productivity improvement compresses marginal cost well below the

pre-shock level and Rotemberg price adjustment passes this through relatively quickly; the

monetary easing (Πt falls) adds a further deflationary component to all nominal sectoral

inflations. Other goods sectors—Agriculture, Construction, and Utilities—also record de-

flation on impact, reflecting the downstream cost relief they receive through IO linkages.

In contrast, several services sectors record mildly positive or near-zero inflation on impact:

while aggregate CPI deflation pushes all nominal sectoral inflations down, the relative price

of services rises with demand-pull pressure from the income expansion, partially or fully

offsetting the aggregate deflation for these sectors. This sign reversal between goods and

services nominal inflation is a signature of a domestically-originated supply shock propa-

gating through the IO network: the cost-relief channel dominates for sectors tightly linked

to the source, while the demand-pull channel dominates for sectors at the periphery of the

Manufacturing subnetwork.

Figure 16 decomposes the aggregate nominal home-price inflation response into expenditure-

weighted goods and services sub-aggregates. Goods nominal home-price inflation falls sub-

stantially on impact, driven by the large deflationary response in Manufacturing itself and

its downstream neighbors; the monetary easing reinforces this, as falling aggregate CPI

(Πt < 1) shifts all nominal sectoral inflations downward. Services nominal home-price infla-

tion is mildly deflationary on impact despite rising relative prices in services: the aggregate

CPI deflation is large enough to more than offset the demand-pull increase in relative service

prices for most services sectors. The expenditure-weighted economy-wide nominal home-

price inflation closely tracks the aggregate CPI response, as it should by construction: the

expenditure-weighted average of 4(∆ log P H

it + log Πt) equals aggregate CPI inflation plus

the consumption-share-weighted average of relative price changes. The remaining wedge be-

tween the home-price aggregate and total CPI reflects the import-price component of the CPI

basket: as the exchange rate depreciates (Q rises), import prices rise in domestic currency,

pushing total CPI slightly above the domestic home-price aggregate. This open-economy

channel operates in the opposite direction relative to the oil shock, where real appreciation

lowered import prices and widened the gap between home-price and CPI inflation in the

60

Figure 15. Sectoral nominal home-price inflation (annualized pp) responses to a 1% Manu-
facturing TFP shock. Inflation is defined as ˆΠH
it−1 (equation (37)), reported as
annualized quarter-over-quarter log changes 4(∆ log P H
it +log Πt). Blue panels: goods sectors
(1–5); red panels: services sectors (6–12).

it = Πt · P H

it /P H

61

other direction.

Figure 16. Expenditure-weighted home-price inflation by sector group following a 1% Man-
ufacturing TFP shock. Aggregate (solid black): economy-wide consumption-share-weighted
average; Goods (dashed blue): within-group average of sectors 1–5; Services (dotted red):
within-group average of sectors 6–12. All in annualized percentage-point deviations from
steady state.

4.4.5 Labor Market Responses

The Manufacturing TFP shock generates an aggregate labor expansion and an increase in

real wages—the mirror image of the oil shock’s labor market contraction. Figure 17 reports

the aggregate employment and wage impulse responses. Employment rises on impact as

Manufacturing directly increases its labor demand to exploit the productivity gain, and

downstream sectors expand production in response to lower intermediate input costs. The

real wage also rises on impact, reflecting the higher marginal product of labor economy-wide:

as Manufacturing becomes more productive, the competitive labor market bids up wages to

equalize the marginal revenue product across all sectors.

Figure 18 displays the sectoral employment responses. The cross-sectoral pattern reveals

the tension between two forces induced by the TFP improvement. Manufacturing employ-

ment itself rises sharply on impact, as the direct productivity gain raises the sector’s labor

demand at given wages. However, with a single competitive wage and quadratic labor ad-

justment costs at the aggregate level, the reallocation of labor toward Manufacturing comes

62

Figure 17. Aggregate labor market responses to a 1% Manufacturing TFP improvement. Left
panel: aggregate employment Nt; right panel: real wage wt. Both in percentage deviations
from steady state.

partly at the expense of other sectors: as the wage rises, sectors that do not benefit directly

from the TFP improvement face higher labor costs without an offsetting productivity gain,

potentially dampening their employment response relative to their output expansion. The

net effect on employment in non-Manufacturing sectors depends on the balance between the

expansionary demand spillover (which raises output and labor demand) and the labor cost

increase (which compresses labor demand at given output). Goods sectors closely linked to

Manufacturing through the IO network—Construction, Utilities, and Agriculture—tend to

expand employment through the cost-relief channel, while services sectors at the periphery

of the Manufacturing subnetwork show more modest employment gains, as the demand-pull

channel is partially offset by higher wage costs.

4.4.6 Output Gap Responses

We now turn to the output gap implied by the Manufacturing TFP shock, which offers a

useful contrast with the oil shock analysis. Figure 19 reports the aggregate output gap and

GDP gap. The positive TFP shock raises the flexible-price output level, but sticky prices

prevent the economy from expanding to its new efficient allocation instantaneously. As a

result, the output gap initially turns negative: actual output undershoots the new, higher

efficient benchmark. Monetary easing accelerates convergence, and the gap closes within

several quarters as prices adjust and the demand expansion catches up to the supply-side

63

Figure 18. Sectoral employment responses to a 1% Manufacturing TFP shock. Manufac-
turing (sector 3) highlighted in red. Blue panels: goods sectors (1–5); red panels: services
sectors (6–12). All series in percentage deviations from steady state.

64

improvement.

Figure 19. Aggregate output gap and GDP gap responses to a 1% Manufacturing TFP
shock. Left panel: output gap (Y gap); right panel: GDP gap. Both in percentage deviations
from steady state.

The sectoral output gaps in Figure 20 reveal a heterogeneous pattern. Manufacturing

itself exhibits a pronounced negative gap—its efficient output expands sharply with the

TFP improvement, but Rotemberg adjustment costs prevent its home price from falling

fast enough to absorb the full demand surge, so actual output falls short of the efficient

level. Downstream sectors with sticky prices also exhibit negative gaps as the deflationary

impulse from cheaper Manufacturing inputs reaches them slowly, distorting relative prices.

Sectors with more flexible prices converge faster, and their output gaps close within a few

quarters. The cross-sectoral dispersion of output gaps underscores that aggregate output

gap measures can be misleading in a multi-sector economy: even when the aggregate gap is

small, individual sectors may face substantial inefficiency wedges driven by the interaction

of heterogeneous price stickiness with IO-transmitted cost changes.

65

Figure 20. Sectoral output gap responses to a 1% Manufacturing TFP shock. Manufacturing
(sector 3) highlighted in red. Blue panels: goods sectors (1–5); red panels: services sectors
(6–12). All series in percentage deviations from steady state.

66

4.5 Comparison: Oil Price Shock vs. Manufacturing TFP Shock

Having characterized the transmission of each shock individually, we now place them side

by side to highlight the distinct and in several respects opposing mechanisms that the NK-

IOSOE model generates. The comparison illuminates how the same input-output network

can act as an amplifier under one type of shock and a redistribution mechanism under

another, and how the direction of monetary policy accommodation shapes the open-economy

offset.

4.5.1 Opposite Aggregate Dynamics

The most striking contrast is at the aggregate level. The oil price shock is stagflationary:

GDP falls by 0.19% on impact while CPI inflation rises by 0.15 annualized percentage points,

a pattern that puts the central bank in the classic dilemma of stabilizing output versus

stabilizing prices. The Manufacturing TFP shock is expansionary-deflationary: GDP rises

by 0.39% and CPI inflation falls by 0.39 annualized percentage points, so the central bank

can ease simultaneously with both output and price objectives pointing in the same direction.

The exchange rate responses are also opposite. Under the oil shock, the central bank

tightens, inducing real appreciation (Q falls by 0.50% on impact): this appreciation is an

open-economy stabilizer that contains the pass-through of oil costs into domestic consumer

prices, but it simultaneously crowds out exports. Under the Manufacturing TFP shock,

monetary easing induces real depreciation (Q rises by 0.19% on impact): the depreciation

partially moderates the domestic deflationary impulse and reinforces export competitiveness.

The trade balance improves under both shocks, but through different mechanisms—import

compression under the oil shock (demand falls due to tightening), and export expansion

under the Manufacturing shock (competitiveness rises due to easing and lower unit costs).

Figure 21 presents the four key aggregate impulse responses under both shocks on a

common set of axes, normalized to the same shock size for comparability.

67

Figure 21. Aggregate impulse responses under the oil price shock (solid orange) and the
Manufacturing TFP shock (dashed blue), normalized to a common shock magnitude. Top
row: GDP and CPI inflation. Bottom row: real exchange rate Q and the policy rate.

68

4.5.2 IO Network as Amplifier vs. Redistributor

Both shocks exploit the input-output network, but the network plays qualitatively different

roles in each case.

Under the oil price shock, the IO network acts as an amplifier. The initial cost increase

in oil-intensive sectors (Transport, Mining, Manufacturing) propagates forward as rising

intermediate input prices for every sector that purchases from them, generating an additional

layer of inflationary pressure beyond the direct cost channel. As shown in Section 4.3, IO

propagation accounts for 60% of the total aggregate marginal cost increase—a factor-of-2.5

amplification over the direct effect alone. Every sector ends up with a higher marginal cost

than it would without IO linkages.

Under the Manufacturing TFP shock, the IO network acts as a redistributor. The pro-

ductivity gain in Manufacturing lowers costs for downstream purchasers of Manufacturing

inputs, easing inflationary pressure in goods sectors and providing a partial offset even for

sectors not directly linked to Manufacturing. At the same time, sectors at the periphery

of the Manufacturing subnetwork—primarily services—face rising prices from demand-pull

forces rather than cost-push forces, as the income and demand effects of the expansion domi-

nate the weak cost-relief channel. The result is a cross-sectoral rebalancing of price pressures

rather than a uniform amplification.

Figure 22 illustrates this contrast directly: the oil shock pushes virtually all sectors to-

ward positive home-price inflation, with a clear ordering by IO centrality and oil intensity;

the Manufacturing TFP shock creates a stark positive-negative split between services (infla-

tionary) and goods (deflationary), with the magnitude of the goods-sector deflation rising

with proximity to Manufacturing in the IO network.

4.5.3 Implications for Monetary Policy

The contrast between the two shocks has direct implications for the conduct of monetary

policy in a network economy.

Under the oil price shock, the Taylor rule induces tightening that dampens inflation but at

the cost of amplifying the output contraction. The IO network magnifies this cost: because

69

Figure 22. Sectoral nominal home-price inflation (annualized pp) under the oil price shock
(orange) and the Manufacturing TFP shock (blue), plotted on the same axes for each sector.
Both series use the definition ˆΠH
it−1 (equation (37)). Blue panels: goods sectors;
red panels: services sectors. Impact response (t = 1) shown as bars; remaining quarters as
shaded bands.

it = Πt·P H

it /P H

70

the inflationary impulse is diffuse—spreading to all sectors through upstream linkages—

the aggregate inflation response is large relative to the originating shock, compelling more

aggressive tightening. The resulting real appreciation provides an open-economy offset to

domestic inflation but exacerbates the output loss.

Under the Manufacturing TFP shock, the policy dilemma disappears: easing simulta-

neously supports output and contains deflation. The IO network redistributes rather than

amplifies the price effects, so the aggregate deflationary pressure that prompts easing is

moderate relative to the size of the supply expansion. The depreciation reinforces the ex-

pansionary channel without risking inflationary overshooting.

More broadly, these results suggest that correctly diagnosing the origin of an inflationary

episode—external cost-push (oil type) versus domestically-generated supply improvement

(TFP type)—matters quantitatively for the appropriate policy response in an economy with

rich IO linkages. A policymaker using a one-sector model would miss the network amplifi-

cation under cost-push shocks and the network redistribution under supply shocks, leading

to systematic misdiagnosis of inflation dynamics and suboptimal policy responses.

5 Conclusion

This paper develops and estimates a New Keynesian model for the Chilean economy with

N = 12 sectors linked by empirical input-output relationships. The model combines hetero-

geneous sectoral price stickiness calibrated from microdata with a three-layer consumption

structure, oil-versus-non-oil disaggregation of intermediate imports, and a debt-elastic small

open economy framework. Structural parameters are estimated by minimizing the distance

between Chilean business cycle moments and their model counterparts using a Simulated

Method of Moments approach.

Our main quantitative finding is that the production network is a first-order amplifier of

shocks in the Chilean economy. For the oil price shock—our leading application—we show

that the input-output multiplier raises the aggregate inflationary impact of a 10% world

oil price increase by a factor of 2.5 relative to the direct cost channel alone: IO network

propagation accounts for 60 percent of the total marginal cost increase, so a model without

71

production linkages would understate the inflationary impact by a factor of 2.5. This am-

plification is heterogeneous across sectors: sectors with moderate own-oil exposure but high

upstream dependence on oil-intensive suppliers experience large second-round cost increases

that a model without IO linkages would miss entirely. The interaction between network

amplification and heterogeneous price stickiness generates rich cross-sectoral dynamics that

are absent from representative-sector models: sectors with stickier prices absorb a larger

fraction of cost shocks through markup compression, dampening their own price response

but amplifying the cost push transmitted downstream.

Several extensions would enrich the framework. Incorporating incomplete exchange-rate

pass-through via Rotemberg pricing in the importable sector, following Romero (2022), would

allow the model to better account for the gradual adjustment of import prices observed in

Chilean data. Allowing for endogenous input-output linkages, as in ?, would let the network

structure respond to relative price changes. On the estimation side, a full Bayesian procedure

with richer data—including wages and trade flows—would improve the identification of labor

market and external sector parameters. We leave these extensions for future work.

References

Ferrante, Francesco, Sebastian Graves, and Matteo Iacoviello, “The inflationary
effects of sectoral reallocation,” Journal of Monetary Economics, 2023, 140, S64–S81.

Gal´ı, Jordi and Tommaso Monacelli, “Monetary Policy and Exchange Rate Volatility
in a Small Open Economy,” The Review of Economic Studies, 07 2005, 72 (3), 707–734.

Pasten, Ernesto, Raphael Schoenle, and Michael Weber, “Sectoral heterogeneity in
nominal price rigidity and the origin of aggregate fluctuations,” Chicago Booth Research
Paper, 2021, (17-25), 2018–54.

Romero, Damian, Domestic linkages and the transmission of commodity price shocks,

Banco Central de Chile, 2022.

, Production Linkages and Nominal Rigidities in a Small Open Economy, Mimeo, 2023.

Silva, Alvaro, Petre Caraiani, Jorge Miranda-Pinto, and Juan Olaya-Agudelo,
“Commodity prices and production networks in small open economies,” Journal of Eco-
nomic Dynamics and Control, 2024, 168, 104968.

72

A Other calculations

The firm solves the following cost-minimization problem (omitting the variety index s):

min
Mt,Lt

PM tMt + PLtLt

s.t. Yt = At

"
α

ϵY −1
ϵY

1
ϵY M
t

+ (1 − α)

ϵY −1
ϵY

1
ϵY L
t

# ϵY

ϵY −1

.

Let M Ct denote the Lagrange multiplier on the production constraint. With this cost-

minimization formulation, M Ct equals the nominal marginal cost (the marginal cost of

producing one additional unit of Yt).

The first-order conditions are

PM t = M Ct

PLt = M Ct

∂Yt
∂Mt
∂Yt
∂Lt

= M Ct AϵY −1

t

(cid:19) 1
ϵY

,

(cid:18)

α

Yt
Mt

= M Ct AϵY −1

t

(cid:18)

(1 − α)

(cid:19) 1
ϵY

.

Yt
Lt

These imply input demands

Mt = α Yt AϵY −1

t

(cid:18) M Ct
PM t

(cid:19)ϵY

,

Lt = (1 − α) Yt AϵY −1

t

(cid:18) M Ct
PLt

(cid:19)ϵY

.

Finally, substituting these expressions back into the production function yields the standard

CES marginal-cost index:

M Ct =

h

1
At

α P 1−ϵY

M t + (1 − α) P 1−ϵY

Lt

i 1
1−ϵY .

B Consumption allocation problems

This appendix collects the intratemporal cost-minimization problems that deliver the de-

mand system in the main text (notation matches the main text).

73

Upper layer: goods vs. services. Given the CES aggregator

Ct =

(cid:16)

ΩSC S
t

σc−1

σc + ΩGC G
t

σc−1
σc

(cid:17) σc

σc−1 ,

the household chooses (C S

t , C G
delivering Ct. The demands are

t ) to minimize nominal expenditure P S

t C S

t + P G

t C G

t subject to

with price index

C g

t = Ωσc
g

!−σc

  P g
t
Pt

Ct,

g ∈ {S, G},

(cid:16)

Pt =

Ωσc

S P S
t

1−σc + Ωσc

G P G
t

1−σc(cid:17) 1

1−σc .

Middle layer: sectors within each category g. For each g ∈ {S, G}, the CES aggre-

gator is

Cost minimization implies

C g

t =

  N
X

i=1

i C g
ωg

it

σI −1
σI

!

σI
σI −1

.

C g

it = (ωg

i )σI

!−σI

  Pit
P g
t

C g
t ,

i = 1, . . . , N,

and

P g

t =

  N
X

i=1

(ωg

i )σI Pit

1−σI

! 1
1−σI

.

Lower layer: home vs. foreign variety for each (i, g). For each i and g, the Armington

aggregator is

(cid:18)

C g

it =

ϱi C Hg
it

σH −1
σH + (1 − ϱi) C F g
it

σH −1
σH

(cid:19) σH
σH −1

.

74

Cost minimization yields

C Hg

it = ϱσH

i

 P H
it
Pit

!−σH

C g
it,

C F g

it = (1 − ϱi)σH

  P F
it
Pit

!−σH

C g
it,

with the corresponding price index

(cid:16)

Pit =

ϱσH
i P H
it

1−σH + (1 − ϱi)σH P F
it

1−σH

(cid:17) 1

1−σH .

C Consumption Structure

D Production Structure

Table 15 describes the nested production technology used by firms in each sector.

75

Table 13. Nested Consumption Structure

Level Aggregator

Description

Parameters

1

Stone-Geary

Total consumption split be-
tween goods and services cat-
egories

ΩG = ¯ω = 0.57

Fixed
shares
expenditure
(non-homothetic preferences)

ΩS = 1 − ¯ω = 0.43

C = ΩG · Cg
Pg

+ ΩS · Cs
Ps

Within-category
across 12 sectors

allocation

γg
i (from data)

Constant expenditure shares
within each category

γs
i (from data)

Cg = Q

i∈Goods C

γg
g,i, Cs = Q
i

i∈Services C

γs
i
s,i

Home vs.
gregation for each sector

foreign variety ag-

ϱi (from data)

Imperfect
substitution be-
tween domestic and imported
goods

σH = 0.999

2

Cobb-Douglas

3

Armington CES

Ci =

Pi =

h
i C(σ−1)/σ
ϱσ
H,i
H,i + (1 − ϱi)σP 1−σ
i P 1−σ
ϱσ

+ (1 − ϱi)σC(σ−1)/σ
F,i
i1/(1−σ)

V

h

iσ/(σ−1)

i (for goods) and γs

Level 1: Households allocate a fixed share of expenditure (57%) to goods and (43%) to services,
independent of relative prices.
Level 2: Within each category, expenditure is allocated across sectors according to Cobb-
Douglas weights γg
i (for services), constructed from Chilean consumption
data. These sum to unity within each category.
Level 3: For each sector i, consumers choose between domestically-produced (CH,i) and im-
ported (CF,i) varieties with home bias parameter ϱi (typically close to 1, indicating strong pref-
erence for domestic goods) and elasticity of substitution σH ≈ 1 (nearly Cobb-Douglas).
Notation: PH,i = domestic price, PV = Q · P ∗
nominal exchange rate.

V = import price in domestic currency, Q =

76

Table 14. Nested Consumption Structure (Compact)

Level Type

Equation

Key Parameters

1

2

3

Stone-Geary

C = ΩG

Cg
Pg

+ ΩS

Cs
Ps

Cobb-Douglas

Cg = Q C

γg
g,i, Cs = Q C
i

γs
i
s,i

ΩG = 0.57, ΩS = 0.43

γg
i , γs

i from data

Armington CES Ci =

h
i C(σ−1)/σ
ϱσ

H,i

+ (1 − ϱi)σC(σ−1)/σ

F,i

iσ/(σ−1)

ϱi from data, σH = 0.999

Level 1: Goods vs.
allocation across 12 sectors. Level 3: Home vs. foreign varieties (Armington aggregation).

services split (fixed expenditure shares). Level 2: Within-category

E Input-Output Accounting Framework

Table 17. Schematic Input-Output Accounting Framework

Intermediate Inputs

Final Demand

S1

S2

S3

S4

S5

S6

Ci

M 1

1t M 2

1t M 3

1t M 4

1t M 5

1t M 6

1t C H
1t

M 1

2t M 2

2t M 3

2t M 4

2t M 5

2t M 6

2t C H
2t

M 1

3t M 2

3t M 3

3t M 4

3t M 5

3t M 6

3t C H
3t

M 1

4t M 2

4t M 3

4t M 4

4t M 5

4t M 6

4t C H
4t

M 1

5t M 2

5t M 3

5t M 4

5t M 5

5t M 6

5t C H
5t

M 1

6t M 2

6t M 3

6t M 4

6t M 5

6t M 6

6t C H
6t

Xi

X1t

X2t

X3t

X4t

X5t

X6t

n
o
i
t
c
u
d
o
r
P

S1

S2

S3

S4

S5

S6

Yi

P

P

P

P

P

P

j M j
j M j
j M j
j M j
j M j
j M j

1t + C H

1t + X1t

2t + C H

2t + X2t

3t + C H

3t + X3t

4t + C H

4t + X4t

5t + C H

5t + X5t

6t + C H

6t + X6t

. Oil
p
m

I

Non-oil

V O
1t

V ¬O
1t

V O
2t

V O
3t

V O
4t

V O
5t

V O
6t

V ¬O
3t

V ¬O
4t

V ¬O
5t

V ¬O
6t

C F
t

V ¬O
2t
t Xt = P

i P X

i PitCit; P X

it Xit; T Bt = P X

Ct = P
Notes: Schematic representation for N = 6 sectors; the calibrated model uses N = 12. M j
it is intermediate
input purchased by sector i from sector j; Vit = V O
it and
Xit are domestic consumption and exports of sector i. The IO coefficient Γi,j gives sector i’s input share
sourced from sector j in the steady state.

is the oil/non-oil import aggregate; C H

i Vi); GDPt = Ct + T Bt

t Xt − P v

it + V ¬O

t + P

t (C F

it

F SMM Estimation

We take the model to Chilean data following a two-step strategy closely related to Ferrante

et al. (2023).

In the first step, we fix a set of parameters externally using steady-state

77

Table 15. Nested Production Structure

Level Aggregator Description

Parameters

1

CES

Gross output produced from
three inputs: domestic materials,
imported inputs, and labor

αmi (from data)

Constant elasticity of substitu-
tion between production factors

αvi (from data)

Labor share is residual: 1−αmi −
αvi

ϵY = 1.484

Yi = exp(Ai)

h
mi M (ϵY −1)/ϵY
α1/ϵY

i

+ α1/ϵY

vi V (ϵY −1)/ϵY

i

+ (1 − αmi − αvi)1/ϵY L(ϵY −1)/ϵY

i

iϵY /(ϵY −1)

2

CES

Domestic materials bundle ag-
gregates
inputs
intermediate
from all 12 sectors

Γi,j (from IO tables)

Input-output linkages determine
sectoral interdependence

ϵm = 0.051

Mi =

hP12

j=1 Γ1/ϵm

i,j

(Mij)(ϵm−1)/ϵm

iϵm/(ϵm−1)

3

Price Index

Materials price index aggregates
prices of all domestic inputs

Γi,j (from IO tables)

Reflects cost of purchasing inter-
mediate inputs

ˆϵm = 0.051

P M

i =

hP12

j=1 Γ1/ϵm

i,j

(P H

j )1−ϵm

i1/(1−ϵm)

Level 1: Each sector produces gross output Yi by combining domestic intermediate inputs
(Mi), imported inputs (Vi), and labor (Li) using a CES technology. The estimated elasticity
ˆϵY = 1.484 > 1 implies production inputs are mild gross substitutes. Shares αmi and αvi vary
by sector (Table 6).
Level 2: The domestic materials bundle Mi aggregates intermediate inputs from all 12 sectors
with weights determined by the input-output matrix Γi,j (Table 8). The estimated elasticity
ˆϵm = 0.051 ≪ 1 implies strong complementarity between IO-linked input varieties, amplifying
network propagation effects.
Level 3: The materials price index P M
reflects the cost of the intermediate input bundle,
weighted by input-output linkages. Sector i purchases inputs from sector j at domestic price
P H
j .
Notation: Ai = sector-specific TFP, Mij = inputs from sector j used by sector i, P H
price of sector j output.
Key difference from consumption: Production has strong network effects through input-
output linkages (Γi,j), while consumption has categorical structure (goods vs. services).

j = domestic

i

78

Table 16. Nested Production Structure (Compact)

Level Type

Equation

Key Parameters

1

2

3

CES

Yi = exp(Ai)

hP

f α1/ϵY

f X (ϵY −1)/ϵY

f

iϵY /(ϵY −1)

αmi, αvi, ˆϵY = 1.484

CES (IO) Mi =

hP12

j=1 Γ1/ϵm

i,j M (ϵm−1)/ϵm

ij

iϵm/(ϵm−1)

Γi,j from data, ˆϵm = 0.051

Price Index P M

i =

hP12

j=1 Γ1/ϵm

i,j

(P H

j )1−ϵm

i1/(1−ϵm)

Network linkages

Level 1: Output from materials (Mi), imports (Vi), labor (Li); estimated elasticity ˆϵY = 1.484.
Level 2: Materials bundle from 12 sectors (IO matrix); estimated elasticity ˆϵm = 0.051 (strong
complementarity). Level 3: Materials cost index. The low IO elasticity ˆϵm ≪ 1 amplifies
shocks through production networks.

relationships and prior literature, as described above. In the second step, we estimate the

remaining structural parameters by minimizing the distance between a set of unconditional

business cycle moments in the data and their model counterparts.

F.1 Shocks

We allow the economy to be hit by three types of shocks. First, a goods-services reallocation

shock, ωt, that alters the relative demand for goods and services. Second, sectoral TFP

shocks, Ai

t, that generate heterogeneous supply-side fluctuations across sectors. Third, a

foreign demand shock, ψit, capturing shifts in external demand for Chilean exports across

sectors, reflecting the SOE nature of the economy.

We assume that each driving term follows an AR(1) process:

ωt+1 = (1 − ρω)¯ω + ρωωt + σωεω
t ,

Ai

t+1 = ρAAi

t + σi

AεA,i

t

,

ψit+1 = (1 − ρψ) ¯ψi + ρψψit + σψεψ,i

t

.

(42)

(43)

(44)

The persistence parameters ρω, ρA, and ρψ and the shock standard deviations σω, σi

A, σψ are

jointly estimated alongside the structural parameters.

79

F.2 Estimated Parameters

We estimate the following structural parameters and shock processes:

• c: the aggregate labor adjustment cost governing hiring frictions across sectors;

• εY : the elasticity of substitution between production inputs (labor, domestic materials,

and imported inputs);

• εm: the elasticity of substitution between sectoral intermediate inputs;

• κv: the price adjustment cost in the importable sector, governing the degree of exchange

rate pass-through to domestic prices;

• ρω, σω: persistence and standard deviation of the goods-services reallocation shock;

• ρA, σA: common persistence and standard deviation of sectoral TFP shocks;

• ρψ, σψ: persistence and standard deviation of the foreign demand shock.

We group these in the vector θ and estimate them by minimizing the distance between a set

of model moments and their data counterparts. Consequently, Tables 3 reports the estimated

values in place of the preliminary calibrated values used for the impulse response analysis in

Section 4.

F.3 Moment Conditions

Rather than targeting changes around a specific episode, we target unconditional second

moments computed from the full sample of quarterly Chilean data (2003Q1–2023Q4). All

series are HP-filtered with smoothing parameter 1600 to extract cyclical components, which

we denote with a hat.

Cross-sectional moments. For each of the N = 12 sectors, we compute the volatility

(standard deviation of the cyclical component) of gross output, the sectoral price defla-

tor, and employment. We denote the resulting cross-sectional vectors by yd, pd, ld ∈ RN

and compute six cross-sectional statistics separately for goods-producing (sectors 1–5) and

80

service-producing (sectors 6–12) groups:

(cid:16)
std

ˆyg(cid:17)

,

(cid:16)
std

ˆpg(cid:17)

,

(cid:16)ˆlg(cid:17)

std

,

(cid:16)
std

ˆys(cid:17)

,

(cid:16)
std

ˆps(cid:17)

,

(cid:16)ˆls(cid:17)

std

,

where std(·) denotes the (output-weighted) cross-sectional average of time-series standard

deviations. We also target the rank correlations between sectoral volatilities in the data and

the model, which discipline how the structure of the IO network and heterogeneous price

stickiness allocate fluctuations across sectors:

1 − ϱ

(cid:16)

(cid:17)
yd, ym(θ)

,

1 − ϱ

(cid:16)

(cid:17)
pd, pm(θ)

,

1 − ϱ

(cid:16)

ld, lm(θ)

(cid:17)

.

Aggregate moments. We target aggregate moments. First, the volatility of aggregate out-
put, std( ˆGDP ). Second, the volatility of aggregate inflation, std(ˆπ). Third, the correlation
between aggregate output and inflation, ϱ( ˆGDP , ˆπ), which captures the model’s ability to

reproduce the slope of the Phillips curve. Fourth, the average goods expenditure share in

consumption, ¯ωG = P G

t C G

t /PtCt, and the services expenditure share ¯ωS = 1 − ¯ωG, which pin

down the steady-state consumption composition. Fifth, the volatility of the real exchange
rate, std( ˆQ), which disciplines the degree of exchange rate pass-through to domestic inflation

through imported inputs. Sixth, the average trade balance-to-GDP ratio, T B/GDP , which

anchors the steady-state external position of the economy.

F.4 Estimator

The estimator solves:

ˆθ = arg min

θ

ψ(θ)⊤W ψ(θ),

(45)

81

where the moment vector ψ(θ) stacks the differences between each data moment and its

model counterpart:





ψ(θ) =

.

(46)









































































ϱ( ˆGDP





















std(ˆyg,d) − std(ˆyg,m(θ))

std(ˆpg,d) − std(ˆpg,m(θ))

std(ˆlg,d) − std(ˆlg,m(θ))

std(ˆys,d) − std(ˆys,m(θ))

std(ˆps,d) − std(ˆps,m(θ))

std(ˆls,d) − std(ˆls,m(θ))

(cid:16)
1 − ϱ

yd, ym(θ)

(cid:17)

1 − ϱ

(cid:16)

(cid:17)
pd, pm(θ)

1 − ϱ

(cid:16)

(cid:17)
ld, lm(θ)

std( ˆGDP

d

) − std( ˆGDP

m

(θ))

std(ˆπd) − std(ˆπm(θ))

, ˆπm(θ))

d

, ˆπd) − ϱ( ˆGDP

m

¯ωG,d − ¯ωG,m(θ)

std( ˆQd) − std( ˆQm(θ))

(T B/GDP )

d

− (T B/GDP )

m

(θ)

82






























































































W is the identity matrix, so all moments receive equal weight. Standard errors are computed

by bootstrap.

F.5 Identification

The parameters are identified through the following logic. The labor adjustment cost c

governs the degree to which reallocation shocks generate persistent relative price movements:

without hiring costs, the reallocation shock would produce no inflation differential between

goods and services sectors, so c is primarily identified by the cross-sectional dispersion of

price volatilities std(ˆpg) and std(ˆps) and their ratio. The production function elasticities εY

and εm govern how shocks propagate through the IO network and affect the relative volatility

of output, prices, and employment across sectors; they are disciplined by the rank-correlation

moments and the within-group dispersion statistics. The shock process parameters (ρω, σω)

are identified by the volatility of the goods expenditure share and its persistence, while

(ρA, σA) are disciplined by aggregate GDP and inflation volatility. The foreign demand shock

parameters (ρψ, σψ) are identified by the volatility of the real exchange rate and the trade

balance, which are informative about the amplitude and persistence of external fluctuations

in a way that domestic shocks alone cannot replicate.

F.6 Data Sources

Cross-sectional output and price data are obtained from the Central Bank of Chile’s quarterly

GDP-by-industry series (Anuario CCNN 2023), which provides gross output and chain-type

price deflators at the 12-sector level. Employment data by sector come from the INE’s

National Employment Survey (ENE). Aggregate inflation corresponds to the Chilean CPI

published by INE. The real exchange rate and trade balance data are from the Central Bank

of Chile’s Balance of Payments statistics.

83

G CES Calibration: Mapping Model Parameters to

Data

This appendix provides a detailed derivation of the mapping between CES (Constant Elas-

ticity of Substitution) function parameters and expenditure shares observed in the data.

Understanding this mapping is crucial for calibrating multi-sector models to match input-

output tables, consumption data, and trade flows.

G.1 General CES Theory and Expenditure Shares

The CES aggregator. Consider a general CES function aggregating n inputs {X1, . . . , Xn}:

Y =

" n
X

i=1

i X (ε−1)/ε
α1/ε

i

#ε/(ε−1)

,

(47)

where αi > 0 are share parameters with Pn

i=1 αi = 1, and ε > 0 is the elasticity of substitution

between inputs.

Cost minimization problem. The producer chooses inputs {X1, . . . , Xn} to minimize

total cost:

min
X1,...,Xn

n
X

i=1

PiXi

subject to Y =

i X (ε−1)/ε
α1/ε

i

#ε/(ε−1)

.

" n
X

i=1

(48)

First-order conditions. Let λ be the Lagrange multiplier (marginal cost). The FOCs

are:

Rearranging:

Pi = λ ·

∂Y
∂Xi

= λ · α1/ε

i

(cid:19)1/ε

.

(cid:18) Y
Xi

Xi
Y

= αi

!ε

.

  λ
Pi

(49)

(50)

84

Price index. Summing over all inputs and using P

i αi = 1:

n
X

1 =

Xi
Y

i=1
" n
X

=

n
X

i=1

αi

  λ
Pi
#1/(1−ε)

!ε

= λε

n
X

i=1

αiP −ε
i

,

λ =

αiP 1−ε
i

.

i=1

This is the unit cost function or price index for output Y .

Expenditure shares. The share of total expenditure on input i is:

si =

=

P

=

PiXi
j PjXj

PiXi
λY
Pi · Y · αi(λ/Pi)ε
λY

= αi

(cid:19)1−ε

.

(cid:18) Pi
λ

Key result: The Cobb-Douglas case. When ε = 1:

si = αi

(exactly, for all prices).

(51)

(52)

(53)

(54)

Implication: For Cobb-Douglas technologies, CES share parameters exactly equal expen-

diture shares, regardless of relative prices. This is why Cobb-Douglas calibration is straight-

forward: simply set αi equal to the observed expenditure share.

Key result: Normalized prices.

If all prices are normalized to unity (Pi = 1 for all i),

then λ = 1 and:

si = αi

(exactly, when Pi = 1).

(55)

Implication: In steady state with all prices normalized to 1, CES parameters equal expen-

diture shares for any ε.

General CES calibration strategy. For ε ̸= 1:

1. If steady-state prices are close to normalization (or similar across inputs), use:

αi ≈ sdata

i

.

85

(56)

2. For large deviations from ε = 1 or substantial price dispersion, invert equation (53):

αi = sdata

i

 P data
i
λdata

!ε−1

,

(57)

where λdata is the observed unit cost.

G.2 Production Function: Labor, Materials, and Imports

Technology. For sector i, output is produced via:

Yi = Ai

h

α1/εY
mi M (εY −1)/εY

i

where:

+ α1/εY

vi V (εY −1)/εY

i

+ (1 − αmi − αvi)1/εY L(εY −1)/εY

i

iεY /(εY −1)

, (58)

• Mi = domestic intermediate inputs,

• Vi = imported intermediate inputs,

• Li = labor,

• αmi = domestic materials share parameter,

• αvi = import share parameter,

• (1 − αmi − αvi) = labor share parameter (residual),

• εY = elasticity of substitution in production.

Cost minimization. The firm solves:

min
Mi,Vi,Li

i Mi + P V Vi + wLi
P M

subject to Yi = Ai[· · · ].

(59)

86

First-order conditions. Let M Ci denote the Lagrange multiplier (marginal cost):

i = M Ci · AεY −1
P M

i

P V = M Ci · AεY −1

i

w = M Ci · AεY −1

i

(cid:18)

(cid:18)

(cid:18)

(cid:19)1/εY

,

(cid:19)1/εY

,

αmi

αvi

Yi
Mi
Yi
Vi

(1 − αmi − αvi)

(cid:19)1/εY

.

Yi
Li

Input demands. Solving for input demands:

Mi = αmiYiAεY −1

i

Vi = αviYiAεY −1

i

!εY

(cid:19)εY

 M Ci
P M
i
(cid:18) M Ci
P V

,

,

Li = (1 − αmi − αvi)YiAεY −1

i

(cid:18) M Ci
w

(cid:19)εY

.

(60)

(61)

(62)

(63)

(64)

(65)

Marginal cost. Substituting back into the production function yields:

M Ci =

1
Ai

h
αmi(P M

i )1−εY + αvi(P V )1−εY + (1 − αmi − αvi)w1−εY

i1/(1−εY )

.

(66)

Expenditure shares. Total nominal cost is:

T Ci = P M

i Mi + P V Vi + wLi = M CiYi.

(67)

The expenditure share on materials is:

smi =

P M
i Mi
M CiYi

P M
i

=

· αmiYiAεY −1

(M Ci/P M

i )εY

i
M CiYi

= αmiAεY −1

i

  P M
i
M Ci

!1−εY

.

(68)

87

Similarly:

svi = αviAεY −1

i

  P V
M Ci

!1−εY

,

sli = (1 − αmi − αvi)AεY −1

i

(cid:18) w

(cid:19)1−εY

M Ci

.

(69)

(70)

Steady-state normalization.

In steady state with Ai = 1 and all prices close to unity:

smi ≈ αmi,

svi ≈ αvi,

sli ≈ 1 − αmi − αvi.

(71)

Calibration from data. From the input-output table and firm-level data, compute:

sdata
mi =

sdata
vi =

sdata
li =

Total domestic intermediate expenditure by sector i
Total costs of sector i
Total import expenditure by sector i
Total costs of sector i
Total labor compensation in sector i
Total costs of sector i

.

,

Then set:

αmi = sdata
mi

, αvi = sdata

vi

.

The labor share is automatically (1 − αmi − αvi) = sdata

li

.

,

(72)

(73)

(74)

(75)

Model implementation.

In main SOE gap.jl, shares are read directly from the cali-

brated data file and assigned as modalpha (intermediate share) and modalphaV (import

share). The labor share is implicit as 1 − αmi − αvi.

G.3 Intermediate Input Aggregator

Technology. Sector i combines inputs from all sectors j = 1, . . . , n into a composite in-

termediate bundle:

where:

Mi =





n
X

j=1



εm/(εm−1)

ij M (εm−1)/εm
β1/εm

ji



,

(76)

88

• Mji = amount of good j used as input by sector i,

• βij = share parameter for input j in sector i’s intermediate bundle,

• εm = elasticity of substitution between intermediate inputs,

• Constraint: Pn

j=1 βij = 1 for each i.

Cost minimization. Sector i chooses {Mji}n

j=1 to minimize:

min
{Mji}

n
X

j=1

P H

j Mji

subject to Mi =



εm/(εm−1)

ij M (εm−1)/εm
β1/εm

ji



,

(77)





n
X

j=1

where P H
j

is the price of sector j’s output.

Input demand. The FOC gives:

Mji = βij

!εm

  P M
i
P H
j

Mi,

where the materials price index is:

P M

i =





n
X

j=1



1/(1−εm)

βij(P H

j )1−εm



.

Expenditure share. The expenditure on input j by sector i is:

Expenditureji = P H

j Mji = βij

!1−εm

  P H
j
P M
i

P M

i Mi.

The expenditure share of input j in sector i’s total intermediate spending is:

sji =

P H
j Mji
P M
i Mi

= βij

!1−εm

.

  P H
j
P M
i

Cobb-Douglas case (εm = 1). When εm = 1:

sji = βij

(exactly).

89

(78)

(79)

(80)

(81)

(82)

Normalized prices.

If P H

j = P M
i

for all j (or all prices equal unity):

sji = βij

(exactly).

(83)

Calibration from Input-Output tables. The IO table reports nominal expenditures

Eji where:

• Row j = supplying sector,

• Column i = purchasing sector,

• Eji = expenditure by sector i on inputs from sector j.

Step 1: Normalize by columns. For each purchasing sector i, compute:

βji =

Eji
k=1 Eki

Pn

=

Sector i purchases from j
Total intermediate spending by sector i

.

(84)

This ensures Pn

j=1 βji = 1 for each i.

Step 2: Matrix notation. Let E be the IO matrix with Eji in row j, column i. Then:

β = E ⊘ c,

(85)

where c = [P

k Ek1, . . . , P

k Ekn] and ⊘ denotes element-wise division (column normalization).

Step 3: Transpose for model code. Some model codes (including ours) define βij

with indices (i, j) where i = purchaser and j = supplier. To match this convention:

modbeta(i, j) = βdata

ji =

Eji
k Eki

P

.

(86)

This requires transposing the normalized IO matrix.

Model implementation.

In the Julia driver main SOE gap.jl, the IO table is loaded and

normalized as:

betaio = readdlm(filename_io, ’,’)

# IO table: (supplier, purchaser)

betax = betaio ./ sum(betaio, dims=1) # Column normalization

90

modbeta = betax’

# Transpose: (purchaser, supplier)

Verification. After calibration, check that:

1. Pn

j=1 modbeta(i, j) = 1 for each sector i (columns in original table sum to 1),

2. In steady state, expenditure shares smodel

ji

≈ βdata

ji when prices are normalized.

G.4 Consumption: Goods vs. Services

Preferences.

In our model, the upper-level consumption aggregator uses Cobb-Douglas

(Stone-Geary) preferences:

log Ct = ωg log C g

t + ωs log C s
t ,

where ωg + ωs = 1. This is equivalent to:

Ct = (C g

t )ωg (C s

t )ωs.

(87)

(88)

Expenditure shares. For Cobb-Douglas preferences, expenditure shares are constant:

pg
t C g
t
t C g
pg
t + ps

t C s
t

= ωg,

ps
t C s
t
t C g
pg
t + ps

t C s
t

= ωs.

This gives demand functions:

t C g
pg

t = ωgPtCt,

t C s
ps

t = ωsPtCt,

where Pt = (pg

t )ωg (ps

t )ωs is the price index.

Calibration. From consumption expenditure data:

ωdata

g =

Total expenditure on goods sectors
Total consumption expenditure

,

s = 1 − ωdata
ωdata

g

.

91

(89)

(90)

(91)

(92)

(93)

(94)

(95)

Model implementation.

In main SOE gap.jl:

ombar_val = 0.57

# Goods share (_g)

# Services share = 1 - ombar_val = 0.43

G.5 Within-Category Consumption: Sectoral Shares

Technology. Within goods and services, consumption follows Cobb-Douglas:

t = Y
C g

(C g

it)γg

i , X

i∈Goods

i∈Goods

t = Y
C s

(C s

it)γs

i , X

γg
i = 1,

γs
i = 1.

Expenditure shares. By Cobb-Douglas properties:

i∈Services

i∈Services

γg
i =

PiC g
it
t C g
pg

t

,

γs
i =

PiC s
it
t C s
ps
t

.

Calibration. From consumption data:

• Classify each sector as goods or services based on expenditure patterns,

• For goods sectors: γg

i = Expenditure on sector i
Total goods expenditure ,

• For services sectors: γs

i = Expenditure on sector i

Total services expenditure.

Model implementation.

In main SOE gap.jl:

goods

= spend_good .> spend_serv

services = spend_serv .> spend_good

gammag

= spend_good ./ sum(spend_good)

gammas

= spend_serv ./ sum(spend_serv)

92

G.6 Home vs. Foreign Varieties (Armington Aggregation)

Technology. For each sector/good i, consumers aggregate home and foreign varieties:

Ci =

h

ϱ1/σH
i

where:

(C H

i )(σH −1)/σH + (1 − ϱi)1/σH (C F

i )(σH −1)/σH

iσH /(σH −1)

,

(96)

• C H

i = consumption of home variety,

• C F

i = consumption of foreign variety,

• ϱi = home bias parameter,

• σH = Armington elasticity.

Demand functions. Cost minimization gives:

C H

i = ϱσH

i

!−σH

 P H
i
Pi

C F

i = (1 − ϱi)σH

  P V
Pi

with price index:

Ci,

!−σH

Ci,

Pi =

h

ϱσH
i

(P H

i )1−σH + (1 − ϱi)σH (P V )1−σH

i1/(1−σH )

.

Expenditure share on home variety. The home consumption share is:

sH
i =

i C H
P H
i
PiCi

= ϱσH
i

!1−σH

.

  P H
i
Pi

(97)

(98)

(99)

(100)

Calibration challenge. Unlike previous cases, ϱi cannot be computed directly from data

because:

1. We observe trade flows, but ϱi depends on unobserved consumer preferences,

2. Relative prices (P H

i /P V ) are endogenous and determined in equilibrium,

93

3. Multiple {ϱi} and price vectors can rationalize the same trade flows.

Joint calibration approach. The model jointly solves for:

• Home bias parameters {ϱi}n

i=1,

• Steady-state prices {P H
i

, Pi, P V , Q},

• Quantities {C H

i , C F

i , Yi, . . .},

subject to matching observed targets:

1. Trade balance / GDP: T B/GDP = targetdata,

2. Home consumption shares by sector: (P H

i C H

i )/(PiCi) = sH,data

i

,

3. Export shares: (P XX)/GDP = targetdata.

Steady-state system. The calibration routine steady ntwsoe calib.m solves a nonlinear

system for:

x = [pH,1, . . . , pH,n, w, Q, C, ϱ1, . . . , ϱn, ωX]⊤,

(101)

with residual equations:

Residuali : Yi − DDi = 0,

i = 1, . . . , n,

Residualn+1 : T B/GDP − targetdata = 0,
n
X

Residualn+2 : N −

Li = 0,

i=1

Residualn+3 : 1 − (pg)ωg (ps)ωs = 0,

Residualn+4:2n+3 : sH,model

i

− sH,data
i

= 0,

i = 1, . . . , n.

(102)

(103)

(104)

(105)

(106)

Implementation. The steady state is solved numerically in the Julia driver main SOE gap.jl

using a nonlinear solver on the full system. The home bias parameters are extracted from

the solution and stored in params nt.modvarrho val.

94

Interpretation. After solving, ϱmodel

i

represents the home bias that, together with equilib-

rium prices, rationalizes observed trade flows and expenditure shares. This is a structural

calibration that respects equilibrium consistency.

G.7 Summary: Calibration Strategy

Table 18 summarizes the mapping from model parameters to data for each CES aggregator.

Table 18. Summary: CES Parameter Calibration

CES Function

Parameter

Elasticity

Calibration Strategy

Production

αmi, αvi

ˆεY = 1.484

Intermediates

βij

ˆεm = 0.051

Direct:
to cost
Set equal
shares from IO table (valid at
normalized SS prices)

Direct: Column-normalize IO
matrix, transpose

Goods vs. Ser-
vices

ωg, ωs

(Cobb-

1
Douglas)

Direct: Set equal to expendi-
ture shares

Within Cate-
gory

γg
i , γs
i

(Cobb-

1
Douglas)

Direct: Normalize sectoral
spending within category

Home vs. For-
eign

ϱi

σH ≈ 1

Jointly solve steady state to
match trade flows

Direct calibration: For elasticities close to 1 (Cobb-Douglas), CES parameters equal expen-
diture shares exactly. Set αi = sdata
Joint calibration: Home bias parameters {ϱi} cannot be observed directly and are solved
jointly with prices to match observed trade patterns and expenditure shares.
Why it works: At normalized steady-state prices (Pi = 1 for all i), the CES expenditure share
equals the CES parameter exactly regardless of the elasticity, so setting αi = sdata
is valid even
for the estimated elasticities (ˆεY = 1.484, ˆεm = 0.051).

.

i

i

G.8 Verification and Diagnostics

After calibration, steady-state cost shares implied by the model are verified to match their

data counterparts: production input shares (smodel

mi ≈ αmi, smodel

vi

≈ αvi), IO expenditure

95

shares (smodel

ji

≈ βij), goods-to-services consumption split (¯ωG,model ≈ 0.57), and trade bal-

ance ratio (T B/GDP model ≈ −0.02). This cross-validation confirms that the CES parame-

terization faithfully represents Chilean expenditure patterns at the calibrated steady state.

96

