The Phillips and Beveridge Curves in a Multi-Sector
Economy

Finn Schüle∗

March 23, 2026

Latest version available here.

Abstract

I develop a New Keynesian model with input-output linkages, search-and-matching

frictions, and sticky prices to study how sector-specific shocks affect output and infla-

tion. Firms draw workers from a common labor pool, creating a novel labor market

propagation channel: higher demand in one sector raises wages and job-finding rates

there, redirects job search, and increases hiring costs elsewhere. Solving the model

nonlinearly shows that sectors’ importance for inflation and monetary policy is state-

dependent—sectors with tighter labor markets raise prices more in response to demand

or supply changes. The resulting Phillips curve steepens as tightness (vacancies over

unemployment) rises, consistent with recent evidence, implying that monetary policy

has larger effects on inflation but smaller effects on output when some sectors are

tight. Calibrated to U.S. data, the model shows that the post-COVID shift toward

goods demand alone explains about 2.9 percentage points of inflation and much of the

observed decline in matching eﬀiciency.

∗Brown University. Email: finn_schuele@brown.edu. I am grateful for the continued guidance and sup-
port of my committee, Gauti Eggertsson, Şebnem Kalemli-Özcan, Stefano Eusepi, and Alexandre Gaillard.
I would also like to thank Fernando Duarte, Amy Handlan, Spencer Kwon, Pascal Michaillat, Siddhi Doshi,
Haoyu Sheng, Matthew DeHaven, Olivia Latus, and other participants at internal lunches and seminars at
Brown for their helpful feedback. This research was conducted using resources and services at the Center
for Computation and Visualization, Brown University. All errors are my own.

1

1 Introduction

In this paper, I introduce search-and-matching frictions into a model with a production
network and sticky prices. The existing production network literature—which models the
economy as an interconnected network of sectors that trade intermediate inputs—shows
that input-output linkages are central to the propagation of sector-specific shocks and their
aggregate effects. Micro shocks—changes in productivity or demand in one sector—create
spillovers to other sectors through prices and quantities.
Intuitively, the macroeconomic
impact of a micro shock on inflation and output depends on the affected sector’s importance
as both a supplier and a user of intermediate inputs. Firms, however, interact not only
through intermediate inputs but also through competition for workers in a shared labor
market. This paper’s main innovation is to incorporate search-and-matching frictions into a
production-network model, allowing me to capture the consequences of these labor market
interactions for aggregate outcomes and their implications for monetary policy.

Incorporating search-and-matching frictions yields three key insights. First, the endoge-
nous reallocation of workers across sectors generates new spillovers through hiring costs. As
labor demand in a sector rises, wages increase, and it becomes easier for workers to find
jobs there. Households respond by searching more in high-wage, high-opportunity sectors
and less in others. This reallocation not only reduces labor supply in other sectors, but also
raises hiring costs across the economy, and alters how sectoral shocks propagate through the
network. In effect, a demand shock in one sector behaves like a negative supply shock in
others, reducing labor supply and raising hiring costs. A key implication is that a sectoral
shock can lower measured matching eﬀiciency—which captures how effectively vacancies and
unemployed workers are converted into new hires—in the economy as a whole, even when
matching eﬀiciency remains constant at the sector level.

Second, search-and-matching frictions generate a nonlinear Phillips curve at the sector
level that steepens with labor market tightness—the ratio of vacancies to unemployed work-
ers. As labor market tightness in a sector rises, it becomes increasingly diﬀicult for firms
in that sector to hire workers. Firms become progressively more supply constrained, forc-
ing them to respond to further changes in supply or demand by adjusting prices rather
than quantities. These dynamics produce a nonlinear Phillips curve—consistent with recent
empirical evidence (Benigno & Eggertsson, 2023; Gitti, 2024)—but through a mechanism
driven by nonlinearities in hiring costs and the elasticity of labor supply rather than by wage
rigidity, a novel mechanism relative to the existing literature.

Third, in contrast to the existing literature modeling the effects of sectoral shocks on
aggregate inflation in a production network, when labor markets feature search-and-matching

2

frictions, a sector’s importance for aggregate inflation becomes state dependent. Sectors with
tighter labor markets exert a stronger influence on overall inflation. As a result, shifts in
demand across sectors can generate endogenous cost-push shocks in the aggregate even when
aggregate demand remains unchanged and sectors are otherwise identical. Rising tightness
in sectors experiencing higher demand leads to larger price increases than the price declines
in sectors facing lower demand. The result is an increase in aggregate inflation.

These findings have important implications for the conduct and effectiveness of mone-
tary policy. A nonlinear Phillips curve implies that additional positive demand or negative
supply shocks have larger effects on inflation when some sectors’ labor markets are tight. Fur-
thermore, because changes in aggregate demand have larger effects on inflation and smaller
effects on output, monetary policy does as well. Consequently, the effects of monetary pol-
icy are also state-dependent: in a tight labor market, policy is a weaker tool for stimulating
output but has greater power to affect inflation. While the main contribution of my paper is
theoretical, this result has clear implications for the post-COVID inflation surge, when peak
year-on-year CPI inflation reached 8.9 percent in June 2022.

The spike in inflation from early 2021 to mid-2022 coincided with a shift in consump-
tion spending from services to goods and an unprecedented surge in labor market tightness
concentrated in the goods sector. This tightness did not manifest in low unemployment (U )
but rather in high vacancies (V ), motivating the focus on V
U as a key indicator of labor mar-
ket conditions. The Federal Reserve kept interest rates low—the first rate hike occurred in
March 2022—despite rising inflation and a rapidly tightening labor market. Pre-pandemic
estimates suggested a flat Phillips curve: changes in unemployment were associated with
only small increases in inflation (Hazell et al., 2022). This implied that the risks of keep-
ing interest rates low and running the economy hot were small—unemployment could fall
substantially without generating large increases in prices.

My model instead suggests that because the goods sector, in particular, was so tight and
therefore capacity constrained, unemployment—which was still around 5.8 percent when in-
flation began picking up in March 2021—did not accurately capture labor market conditions.
Tightness in the goods sector, on the other hand, was already 45 percent higher than average
in March 2021, and continued rising to nearly 150 percent above average. Keeping the pol-
icy rate low in this environment had a large effect on inflation but only a limited impact on
unemployment—the opposite of what the Fed hoped to achieve. Had they tightened earlier,
they could likely have reduced inflation substantially without causing much unemployment.
I begin by developing a framework that combines three well-established components from
the literature: a production network, sticky prices, and search-and-matching in the labor
market. I then use this framework to think about the effect of sectoral shocks on inflation

3

and monetary policy, one of many possible applications. The production network allows
me to take a sectoral approach to inflation—following Baqaee and Farhi (2022), La’O and
Tahbaz-Salehi (2022), and Rubbo (2023)—by tracing how shocks propagate downstream
to users of a sector’s output and upstream to its input suppliers. Households optimally
send more workers to sectors with higher wages and a higher probability of finding a job,
triggering cascades in hiring costs, and creating a labor reallocation channel absent from
previous models.

The labor reallocation channel is reminiscent of a classic literature debating whether sec-
toral reallocation caused by sector-specific shocks can generate persistent frictional unem-
ployment (Abraham & Katz, 1986; Lilien, 1982). In my framework, frictional unemployment
exists even in the absence of sectoral shocks; it instead arises from search-and-matching fric-
tions as in the Diamond-Mortensen-Pissarides model (Diamond, 1982b; Mortensen, 1982a;
Pissarides, 1984). However, sectoral shocks do cause costly reallocation of labor across sec-
tors, altering the propagation of shocks to marginal costs and, ultimately, affecting measured
aggregate matching eﬀiciency.1

Unlike most prior work on inflation and monetary policy in production networks, I
solve the model nonlinearly by simulating perfect-foresight impulse responses to unexpected
shocks.2 Using an illustrative calibration with identical sectors designed to resemble an av-
erage U.S. sector, I show that nonlinearities are quantitatively important for both inflation
and the effectiveness of monetary policy.

For instance, when labor market tightness is roughly 125% above steady state—an in-
crease from 0.65 to 1.47 in the goods sector—a policy-rate cut has nearly ten times the effect
on inflation within that sector as under steady-state tightness. For comparison, at peak,
tightness in the goods sector reached to 1.60 in December 2022. An 8% shift in relative de-
mand from services to goods—similar in magnitude to the post-pandemic reallocation—raises
aggregate inflation by 7 percentage points, reduces measured aggregate matching eﬀiciency
by 25%, and amplifies the inflationary impact of a policy-rate cut by roughly 150%, even
though aggregate demand remains unchanged.

A linearized version would miss the steepening of the Phillips curve at high levels of tight-
ness, thereby underestimating the effects of aggregate-demand stimulus on inflation when
labor markets are tight and failing to capture the inflationary effects of demand reallocation.
Furthermore, the nonlinear solution shows that search-and-matching frictions not only alter

1Aggregate matching eﬀiciency measures the effectiveness of the matching process; it is effectively a resid-
ual capturing the wedge between the observed number of hires and the number predicted by the aggregate
numbers of unemployed workers and vacancies in the economy.

2Although the perfect-foresight assumption abstracts from risk, it allows me to isolate the nonlinearities

arising from labor market frictions.

4

which sectors are most influential for inflation through the labor market spillover channel—
that would persist even in a linearized model—but also make a sector’s importance state
dependent.

After establishing the key theoretical mechanisms in the illustrative calibration, I apply
an empirical calibration—accounting for differences in size, labor share, and other sectoral
characteristics—to the recent inflation surge, focusing on the distinction between goods
and services.
I calibrate the production network parameters to U.S. input-output tables
from the Bureau of Economic Analysis (BEA) and estimate the remaining labor market
parameters in a linearized version of the model using Bayesian methods and data from 2000-
2019—a period when labor market tightness remained below the exceptionally high post-
pandemic levels, making the linearized model a reasonable approximation. The persistent
rise in the consumption share of goods relative to services after the pandemic increased
inflation by about 2.9 percentage points and explains nearly the entire decline in aggregate
matching eﬀiciency.3 However, the sectoral demand shock alone cannot explain the entire
rise in aggregate inflation and generates a counterfactual decline in labor market tightness in
services. Incorporating shocks to the separation rate and aggregate demand yields a much
closer fit to the data, and only a modest amount of monetary stimulus leads to large effects
on inflation.

My paper relates to three broad strands of the literature. First, it builds a production
network model in the tradition of Acemoglu et al. (2012), Baqaee and Farhi (2019), Baqaee
and Rubbo (2022), V. M. Carvalho and Tahbaz-Salehi (2019), and Jones (2011). Much of this
literature focuses on the effects of sector specific productivity shocks. According to Hulten’s
theorem—a fundamental benchmark in the literature—the aggregate effect of a productivity
shock in a sector is proportional to its Domar weight, the sector’s income over GDP (Hulten,
1978). Like other ineﬀiciencies explored in the literature—for instance markups and rigid
prices—labor market frictions lead to violations of Hulten’s theorem, altering how shocks
propagate and which sectors are most important for aggregate outcomes. The literature
typically accounts for these violations through a model-informed reweighting of the input-
output linkages, leading to altered—but fixed—quasi-Domar weights.

My contribution is twofold. First, I show how incorporating search-and-matching frictions
in the labor market alters the propagation of sectoral shocks, extending Schüle and Sheng
(2024)4. Relative to Schüle and Sheng (2024), I incorporate nominal rigidities and dynamics
to analyze inflation and monetary policy, and I solve the model nonlinearly. The presence

3The difference in magnitude between the 2.9 percentage point number and the 7 percentage points in the
illustrative example comes from the fact that the goods sector is smaller than the service sector, accounting
for only about 33% of consumption, and has a smaller labor share.

4This research agenda has grown out of collaboration with Haoyu Sheng, a fellow graduate student.

5

of labor market frictions implies that a sector’s importance for aggregate inflation depends
on how shocks pass through to sectoral tightness, how sectoral tightness affects hiring costs,
and how those costs transmit to inflation—necessitating a reweighting of sectoral importance
even in the linearized case. Furthermore, I show that sectoral importance for output and
inflation is state dependent. When labor markets feature search-and-matching frictions,
sectors with tighter labor markets experience larger increases in hiring costs for the same
increase in tightness, and pass on more of any change in demand or supply to inflation.
Therefore, even the static reweighting is not suﬀicient to capture a sector’s importance for
inflation: accurately capturing which sectors are likely to drive inflation requires a dynamic
and state-dependent approach.

Second, my paper relates to the extensive literature that uses multi-sector models to
study the post-pandemic economy (Amiti et al., 2023; Baqaee & Farhi, 2022; Boehm &
Pandalai-Nayar, 2022; Comin et al., 2023; Di Giovanni et al., 2023; Ferrante et al., 2023;
Guerrieri et al., 2022). This literature builds on both the production network literature
highlighted above and an existing literature that developes multi-sector or multi-region New
Keynesian models, for instance McLeay and Tenreyro (2020) and Hazell et al. (2022). Guer-
rieri et al. (2022) show that sectoral supply shocks can trigger demand contractions in other
sectors, pushing aggregate demand below potential—a phenomenon they term a “Keynesian
supply shock”. In my model, labor market spillovers reverse this result: an increase in de-
mand in one sector reduces labor supply in others. As a result, even if aggregate demand
remains unchanged, shifts in relative demand across sectors reduce potential output through
labor reallocation, inducing a positive cost-push shock. These results are complementary:
both show that changes in demand or supply in one sector lead to spillovers in others—
demand shocks are no longer purely demand shocks, and supply shocks are no longer purely
supply shocks. The mechanisms differ, however: Guerrieri et al. (2022) emphasize demand
spillovers through income and substitution effects, whereas I focus on labor market spillovers
driven by search-and-matching frictions, which are especially important when shocks lead to
substantial heterogeneity in tightness across sectors, as during the pandemic recovery.

Like Amiti et al. (2023), Comin et al. (2023), Di Giovanni et al. (2023), and Rubbo (2024),
I examine how sectoral shocks affect inflation. My contribution is to show that search-and-
matching frictions generate a nonlinear Phillips curve at the sector level and state-dependent
effects of sectoral shocks. These frictions also amplify the inflationary effects of changes in
productivity or demand when labor markets are tight. Given the unusually tight labor
markets of the post-pandemic period, my mechanism implies that even modest changes in
sectoral supply or aggregate demand can generate substantial inflationary pressure. Conse-
quently, the risk of keeping interest rates low for too long is greater than standard models

6

suggest.

Within this literature, my paper is most closely related to Ferrante et al. (2023) and
Boehm and Pandalai-Nayar (2022). Relative to Ferrante et al. (2023) I microfound hiring
costs using standard search-and-matching, which naturally generates larger hiring costs in
tighter labor markets and makes hiring costs depend not only on a sector’s own hiring deci-
sions but, through labor reallocation, also on the hiring decisions of other sectors. Relative to
Boehm and Pandalai-Nayar (2022), I show that search-and-matching frictions can effectively
induce capacity constraints at the sector level, amplifying the effects of supply and demand
disturbances on inflation and altering the effectiveness of monetary policy.

Third, I contribute to recent work exploring the links between the inflation surge and
labor market conditions.
I build on Blanchard and Galí (2010), who incorporate search-
and-matching into an otherwise standard New Keynesian model, by adding multiple sectors
that interact through a production network. This extension allows me to analyze the effects
of sectoral shocks and to highlight important cross-sector spillovers that shape aggregate
inflation dynamics.

Benigno and Eggertsson (2023) and Michaillat and Saez (2024) propose mechanisms
that generate a nonlinear Phillips curve in labor market tightness—an empirically plausible
feature (Gitti, 2024)—by assuming downwardly rigid nominal wages. By solving the model
nonlinearly, I show that a smooth steepening of the Phillips curve arises from incorporating
search-and-matching into a standard New Keynesian model with Nash bargained wages.
The nonlinearity emerges from the convexity of hiring costs—it becomes progressively more
diﬀicult to fill vacancies as tightness rises. By focusing on a multi-sector economy, I also
show how this nonlinearity induces state-dependent sectoral importance and endogenous
cost-push shocks when demand shifts across sectors.

Furthermore, I show that sectoral labor reallocation and heterogeneity in tightness across
sectors can lead to endogenous declines in aggregate matching eﬀiciency, shifting the Bev-
eridge curve outward, as observed during the pandemic recovery. Workers moving across
sectors can be viewed as a form of job-to-job transition, which Bagga et al. (2025) show helps
explain the outward shift in the Beveridge curve. Allowing for endogenous separations in
response to wage differentials across sectors would strengthen the power of the reallocation
channel I highlight. Similarly, Afrouzi et al. (2024) propose a mechanism linking inflation to
separations, whereby workers leave their jobs when inflation is high to recoup losses in real
earnings. Again, I view my results as complementary: an increase in demand in one sector,
which raises wages and prices there, would induce higher separations and greater labor flows
across firms, further amplifying the labor spillover channel I highlight.

The paper proceeds as follows. Section 2 introduces the model and outlines the equi-

7

librium conditions. Section 3 presents key predictions. Section 4 applies the model to
post-pandemic inflation and labor market dynamics. Section 5 concludes.

2 A Model with Production Networks and labor mar-

ket Frictions

In this section, I build a New Keynesian model with labor market frictions and production
linkages. The model has three key components: (1) sticky prices, (2) a production network,
and (3) search-and-matching frictions in the labor market. Sticky prices generate a New
Keynesian Phillips curve, linking inflation to firms’ marginal costs and expected future in-
flation. The production network—firms use intermediate inputs from other sectors—links
marginal costs in one sector to prices in other sectors. Finally, search-and-matching frictions
capture the costly hiring process, linking firms’ marginal costs to labor market conditions.
Together, these three components generate a sector-level Phillips curve that captures how
inflation spreads across sectors through both prices and labor market conditions.

2.1 The Labor Market

I begin with the labor market, modeled using a search-and-matching framework in which
firms post vacancies and workers search for jobs (Diamond, 1982a, 1982b; Mortensen, 1982a,
1982b; Pissarides, 1984, 1985). The search-and-matching process generates a probability of
a worker finding a job and of a firm filling a vacancy that depend on sectoral labor market
tightness—the ratio of vacancy postings to unattached job seekers in sector i. As a result,
tightness directly affects firms’ hiring costs, and therefore their marginal costs.

At the start of each discrete time period t, firms in sector i ∈ {1, · · · , J} post vacancies
Vi,t and hire Hi,t workers from a pool of searchers Ui,t.5 A hire occurs whenever a searching
worker matches with a hiring firm. As is standard in the search-and-matching literature, I
represent this matching process using a matching function mi,t(Ui,t, Vi,t).6 The total number
of hires in sector i is

Hi,t = mi,t(Ui,t, Vi,t).

(1)

5While I denote total searchers by Ui,t for consistency with the notation commonly used in the search-
and-matching literature, I think of Ui,t as capturing the number of unattached workers searching for a job
in sector i.

6By allowing the function mi,t to be sector-specific and time varying, I can allow for exogenous variation

in the eﬀiciency of the matching process across sectors and over time.

8

The matching function exhibits constant returns to scale, mi,t(λUi,t, λVi,t) = λHi,t, is increas-
> 0, and concave, ∂2mi,t
ing in both arguments, ∂mi,t
< 0. Lastly, the matching
∂U 2
∂Ui,t
i,t
function satisfies mi,t(Ui,t, Vi,t) ≤ min{Ui,t, Vi,t}. This condition ensures that the numbers of
vacancies and unattached workers at the end of the period t are both non-negative. Intu-
itively, it implies firms cannot conjure workers out of thin air by posting more vacancies.

, ∂2mi,t
∂V 2
i,t

, ∂mi,t
∂Vi,t

Because the matching function is constant returns to scale, I can conveniently express
the probability of a firm filling a vacancy (or the vacancy-filling rate) in sector i in terms of
sector-specific labor market tightness θi,t = Vi,t
Ui,t

:

qi,t =

Hi,t
Vi,t

= mi,t

(cid:0)

(cid:1)

.

θ−1
i,t , 1

(2)

Since mi,t is increasing in both arguments, the vacancy-filling rate falls as tightness rises: it
gets harder for firms to hire as tightness rises. Similarly, the probability that a worker finds
a job in sector i (or the job-finding rate), conditional on searching in that sector, is

fi,t =

Hi,t
Ui,t

= mi,t(1, θi,t).

(3)

The job-finding rate is increasing in tightness:
it gets easier for workers to find a job as
tightness rises. Since mi,t(Ui,t, Vi,t) ≤ min{Ui,t, Vi,t}, both the vacancy-filling rate and the
job-finding rate are bounded between 0 and 1: qi,t, fi,t ∈ [0, 1].

I assume that a vacancy requires ri,t hours of labor to post and maintain—firms pay a
recruiting cost in terms of labor, representing time spent on posting vacancies, reviewing
applications, interviewing candidates, and other recruiting activities. As tightness rises,
firms must post more vacancies, employ more recruiters, and incur higher recruiting costs
to achieve a given level of hiring. This creates a link between labor market tightness and
marginal costs. Finally, I assume that an exogenous fraction si,t of workers separate from
their jobs at the beginning of each period.

2.2 Firms

The production structure in each sector incorporates the three key model ingredients: (1)
sticky prices, (2) a production network, and (3) search-and-matching frictions. To allow
for sticky prices, I assume that each sector contains a continuum of monopolistically com-
petitive firms that use labor and intermediate goods to produce output and set prices sub-
ject to adjustment costs. Price adjustment costs generate a sector-specific New Keynesian
Phillips curve linking inflation to marginal costs and expected future inflation. As in Baqaee

9

(2018),Baqaee and Farhi (2019), and Baqaee and Rubbo (2022), firms use intermediate
inputs from other sectors, generating a production network. Finally, firms make hiring de-
cisions subject to the search-and-matching structure described above, linking marginal costs
to labor market conditions.

Each sector contains a unit continuum of monopolistically competitive intermediate goods
producers indexed by z, and a representative firm that aggregates the intermediate goods to
produce sectoral output. The representative firm in each sector produces the sectoral output
Yi,t from the outputs of each of the intermediate firms Yi,t(z) using CES technology.

Yi,t =

(cid:18)Z

1

0

(cid:19) ϵ
ϵ−1

Yi,t(z)

ϵ−1
ϵ

(4)

where ϵ is the elasticity of substitution between the outputs of the different intermediate
goods producers in sectoral output. This representative firm chooses how much of each firm
z’s output to use in production to minimize its costs, resulting in downward-sloping demand
functions for each z:

Yi,t(z) =

(cid:19)−ϵ

(cid:18)

Pi,t(z)
Pi,t

Yi,t

(5)

where Pi,t(z) is the price of firm z’s output.

As in the standard production networks setup, the monopolistically competitive inter-
mediate goods producers in each sector produce their output Yi,t(z) using labor Ni,t(z) and
intermediate inputs Xij,t(z).

Yi,t(z) = Ai,t

1
ϵy

ω

in Ni,t(z)

ϵy −1
ϵy +

JX

j=1

1
ϵy

ω

ij Xij,t(z)

! ϵy
ϵy −1

ϵy −1
ϵy

(6)

where Ai,t is a sector-specific productivity shock, which captures, for instance, supply chain
disruptions that disproportionately affect one sector. ωin is the share of labor and ωij is the
share of sector j’s output in sector i’s production process. ϵy is the elasticity of substitution
between labor and intermediate inputs.

Each firm z chooses how many vacancies Vi,t(z) to post, taking the total number of
1
vacancy postings in sector i, Vi,t =
0 Vi,t(z)dz, and therefore sector level tightness θi,t, as
given. The number of newly hired workers at firms z, Hi,t(z) depends on the vacancy-filling
rate qi,t, implying the following law of motion for employment at firm z:

R

Ni,t(z) + N r

i,t(z) = qi,tVi,t(z) + (1 − si,t)

(cid:0)

Ni,t−1(z) + N r

i,t−1(z)

(cid:1)

.

(7)

10

Where Ni,t(z) is labor hours engaged in production and N r
spent recruiting.

i,t(z) = ri,tVi,t(z) is labor hours

Finally, each firm z chooses its price Pi,t(z) subject to Rotemberg adjustment costs.7,
the number of vacancy postings Vi,t(z), and how much of each intermediate input to use in
production {Xij,t(z)}J
j=1 in order to maximize its dividend payments8,

max

{Pi,t+s(z),Vi,t+s(z),{Xij,t+s(z)}J

j=1

(cid:20)

∞X

Et

}∞
s=0

s=0

SDFt|t+s

(cid:21)

Di,t+s(z)
Pt+s

(8)

subject to Eqs. (5), (6), and (7). Where

Di,t(z)
Pt

=

−

Pi,t(z)
Pt
JX

(cid:19)−ϵ

(cid:18)

Pi,t(z)
Pi,t

Yi,t − Wi,t
Pt
(cid:18)

Pj,t
Pt

Xij,t(z) − ψp,i
2

Pi,t(z)
ΠPi,t−1(z)

j=1

(cid:2)

Ni,t(z) + N r

i,t(z)

(cid:3)

(cid:19)

2

− 1

Yi,t.

Since all firms in a given sector face an identical problem, they all make identical input
choices and pricing decisions. I therefore drop the z indexation from the following optimality
conditions. I allow the parameter governing the size of price adjustment costs, ψp,i, to vary
at the sector level (Ferrante et al., 2023; Pasten et al., 2020).

Firms post vacancies up to the point where the marginal cost of an additional vacancy

equals the marginal benefit of an additional vacancy

ri,t = µi,t(qi,t − ri,t) + Et

(cid:2)

Wi,t
Pt

SDFt|t+1Πt+1(1 − si,t+1)µi,t+1ri,t

(cid:3)

(9)

where µi,t is the marginal value of an additional employee to the firm. The marginal value
of an additional employee equals the marginal product of labor today, net of the wage, plus
the continuation value of the employee.

µi,t = λi,tβ

1
ϵy

in A

ϵy −1
ϵy
i,t N

− 1
ϵy
i,t Y

1
ϵy
i,t

− Wi,t
Pt

(cid:2)

+ Et

SDFt|t+1Πt+1(1 − si,t+1)µi,t+1

(cid:3)

(10)

where λi,t is the real marginal cost. Each firm chooses intermediate inputs so that good j’s

7I choose Rotemberg adjustment costs because they do not lead to price dispersion across firms, and I

therefore do not need to track the distribution of prices, output, and therefore employment across firms.
8By choosing Vi,t(z) the firm is also choosing Ni,t(z) through the law of motion for employment.

11

share in total costs satisfies

Pj,t
Xij,t
Pt
λi,tYi,t

(cid:18)

= (βixωij)

1
ϵy

Ai,t

(cid:19) ϵy −1
ϵy

Xij,t
Yi,t

(11)

Each firm’s optimal price setting implies a sector-specific Phillips curve as in Ferrante et al.
(2023), La’O and Tahbaz-Salehi (2022), and Rubbo (2023, 2024).
(cid:18)

(cid:18)

(cid:19)

(cid:19)

(cid:18)

(cid:19)

(cid:20)

(cid:21)

Πi,t
Π

− 1

Πi,t
Π

=

ϵ
ϕp,i

M Ci,t
Pt

− ϵ − 1
ϵ

Pi,t
Pt

+ Et

SDFt+1|tΠt+1

Πi,t+1
Π

− 1

Πi,t+1
Π

Yi,t+1
Yi,t
(12)

where SDFt+1|t = β Uc,t+1
Uc,t

t+1 is the household’s stochastic discount factor.
Finally, each sector’s output must satisfy the market clearing condition

Π−1

JX

Yi,t = Ci,t +

Xji,t

j=1

(13)

2.3 Households

A unit mass of identical households indexed by z supply labor to firms, consume goods
produced by each of the J sectors, and save in the form of bonds.

Households have CES preferences for consumption across the sectors. The aggregate

consumption good Ct is therefore:

Ct =

JX

i=1

1
ϵd
i,t C

α

ϵd

−1
ϵd

i,t

! ϵd
−1
ϵd

(14)

where ϵd is the elasticity of substitution across sectors, and αi,t is a time-varying exogenous
preference shifter that captures households’ changing preferences for consumption in sector i.
Note that Ci,t is final household consumption, not total output, because some of the output
from each sector is used as an intermediate input in production.

These preferences imply downward-sloping consumption demand functions:

Ci,t = αi,t

(cid:19)−ϵd

(cid:18)

Pi,t
Pt

Ct

(15)

12

where the aggregate consumer price index is

! 1

1−ϵd

Pt =

JX

i=1

αi,tP 1−ϵd

i,t

(16)

Each household maximizes lifetime utility subject to their budget constraint and the law

of motion for employment in each sector, which depends on the job-finding rates.

max

{Ct+s(z),{Li,t+s(z)}J

i=1}∞

s=0

Et

∞X

s=0

(cid:18)

n

βsUt+s

Ct+s(z),

˜Li,t+s(z), Li,t+s(z)

(cid:19)

o

J

i=1

s.t. Pt+sCt+s(z) + Bt+s(z) = (1 + it+s)Bt+s−1(z) +

JX

i=1

Wi,t+s

˜Li,t+s(z) + Tt+s(z)

˜Li,t+s(z) = (1 − fi,t+s)(1 − si,t+s) ˜Li,t+s−1(z) + fi,t+sLi,t+s(z)
JX

and

i=1

Li,t+s ≤ Lt+s

(17)

where Ct(z) is the aggregate consumption bundle, Li,t(z) is labor supplied to sector i (includ-
ing employed and searching workers), ˜Li,t(z) is employment in sector i, Bt(z) is a one-period
risk-free bond, it is the nominal interest rate, and Tt+s(z) are lump-sum transfers and div-
idend payments. The model features a form of directed search: households choose how to
i=1(1−si,t) ˜Li,t−1(z), across the sectors
divide their endowment of unattached workers, Lt −
based on both the real wage and the likelihood of finding a job in each sector.

P

J

The households’ period utility function is:

Ut+s(·) = Zt+s

"

Ct+s(z)1−σ
1 − σ

−

JX

i=1

χi,t+s

˜Li,t+s(z)1+φ
1 + φ

+

ψL,i
2

(cid:18)

Li,t+s/Lt+s
Li,t+s−1/Lt+s−1

(cid:19)2

!#

− 1

Li,t+s

.

(18)

(cid:16)

(cid:17)

2

Li,t+s/Lt+s
Li,t+s−1/Lt+s−1

The final term in the households’ utility function, ψL,i
2
cost of adjusting the fraction of workers who search in a given sector i. This adjustment
cost captures frictions in relocating across sectors because of moving costs, retraining costs,
or other frictions. These types of frictions are empirically important in determining where
workers search for jobs (Humlum, 2021), and within the model will be important for deter-
mining the labor supply response to sectoral shocks. ψL,i governs the size of reallocation
costs while φ governs how responsive labor supply is to changes in real wages.

− 1

Li,t+s, is a quadratic

The household’s first-order conditions imply the standard Euler equation for consump-

13

tion:

(cid:20)

C −σ

t = β(1 + it)Et

C −σ
t+1

(cid:21)

Zt+1
Zt

Π−1
t+1

(19)

The household optimal decisions for how many workers to send in search of a job in each
sector equalizes the expected return from searching in each sector, subject to the labor
adjustment costs.

δt = Ξi,tfi,t − ψL,iZt

"

(cid:18)

+ βψL,iEt

Zt+1

 (cid:20)

Li,t/Lt
Li,t−1/Lt−1

(cid:21)

− 1

+

(cid:20)

1
2

(cid:19)

!

(cid:21)2

Li,t/Lt
Li,t−1/Lt−1

#

− 1

Li,t+1/Lt+1
Li,t/Lt

− 1

L2

i,t+1/Lt+1
L2
i,t/Lt

(20)

where Ξi,t is the value to the household of an additional employed worker in sector i , fi,t is
the job-finding rate in sector i, and δt is the marginal value of an additional unit of aggregate
labor to the household. Absent labor adjustment costs, the household adjusts search across
sectors up to the point where fi,tΞi,t = fj,tΞj,t. fi,tΞi,t is the expected return from searching
in sector i. It is the probability of finding a job in sector i times the value to the household
of an additional employee in sector i.

The value of an additional employed worker in sector i is the value of being employed
today—how much the employed worker earns, converted into utils, net of the utility cost of
an additional employee—plus the continuation value of being employed tomorrow.

Ξi,t = Zt

(cid:18)

C −σ
t

Wi,t
Pt

(cid:19)

− χi,t

˜Lφ
i,t

+ βEt [(1 − fi,t+1)(1 − si,t+1)Ξi,t+1]

(21)

Finally, the number of unattached workers at the beginning of the period is given by

Ui,t = Li,t − (1 − si,t) ˜Li,t−1.

2.4 Wages

As in standard search-and-matching models, when workers are matched with firms they
face a bilateral monopoly. Wages are therefore not determined by the model’s equilibrium
conditions,9 and instead require an assumption about the wage-setting process.

I use Nash bargaining between the household and the firm as the basis for my wage
i,t denote the firm surplus from a match in sector i and let S h

determination process. Let S f

i,t

9The household and firm equilibrium conditions define a range of acceptable wages for a given match.

14

denote the household surplus from a match in sector i. The Nash bargaining solution implies

i,t = κS f
S h

i,t

(22)

where κ is determined by the relative bargaining power of the household.10

Because hiring an additional worker requires

vacancy postings at cost ri,tWi,t, firms
Wi,t. By free entry into vacancy postings, the firm surplus

1
qi,t

can replace any worker at a cost ri,t
qi,t
from an existing worker is:

S f

i,t =

ri,t
qi,t

Wi,t

(23)

From the household’s first order conditions, the value of an additional employee to the
household, measured in utils, is given by Ξi,t in Eq.(21). The surplus in nominal terms is
therefore Ξi,t
Pt. The Nash bargaining solution then implies that the Nash-bargained wage
ZtC
is

−σ
t

Wi,t =

qi,t
qi,t − κri,t

(cid:18)

(cid:20)

χi,t

˜Lφ
i,tC σ

t Pt − κEt

SDFt+1|t(1 − fi,t+1)(1 − si,t+1)

(cid:21)(cid:19)

Wi,t+1

(24)

ri,t+1
qi,t+1

As is well understood since Shimer (2005) and R. E. Hall (2005), fully flexible Nash bar-
gaining leads to counterfactually large wage fluctuations and small employment fluctuations
over the business cycle. I therefore allow for a degree of wage rigidity by assuming that the
final wage is a weighted average of this Nash bargaining solution and last period’s wage and
I allow for exogenous shocks to the wage, W s

i,t, so that the overall wage is given by

Wi,t =

(cid:18)

(cid:20)

qi,t
qi,t − κri,t

(cid:20)

χi,t

˜Lφ
i,tC σ

t Pt − κEt

SDFt+1|t(1 − fi,t+1)(1 − si,t+1)

(cid:21)(cid:19)(cid:21)

1−ρw

ri,t+1
qi,t+1

Wi,t+1

(25)

× (Wi,t−1)ρw W s
i,t

This specification nests fully flexible Nash bargained wages (ρw = 0) and fully rigid wages
(ρw = 1). This is a reduced-form way to capture that wages may not adjust fully within a
month to changes in economic conditions, and as Gertler and Trigari (2009) demonstrate,
partial wage rigidity allows the search-and-matching model to generate more realistic fluc-
tuations in wages and unemployment.

10If, for instance, the firm and household split the total surplus equally, then κ = 1.

15

2.5 Monetary Policy and Shock Processes

I close the model by assuming that the central bank sets the nominal interest rate according
to a Taylor rule in inflation and the output gap.

(1 + it) = (1 − it−1)ρi

(cid:21)

ϕπ

(cid:20)

 (cid:20)

Πagg
t
Π

Y agg
t
Y agg,f lex
t

!

(cid:21)

ϕy

1−ρi

εm
t

(26)

where εm
t

is a normally distributed monetary policy shock.

Finally, I assume that the natural logarithms of Zt, χi,t, αi,t, si,t, ri,t, W s

i,t, follow AR(1)

processes with normally distributed innovations.

log Zt = ρz log Zt−1 + εz
t ,
log αi,t = ρα log αi,t−1 + εα
i,t,
log ri,t = ρr log ri,t−1 + εr
i,t,

log χi,t = ρχ log χi,t−1 + εχ
i,t
log si,t = ρS log si,t−1 + εs
i,t
i,t−1 + εW s

i,t = ρw log W s

i,t

log W s

(27)

(28)

(29)

I summarize the full set of equilibrium conditions in Appendix A.11

2.6 Equilibrium Definition

An equilibrium is a set of sequences for endogenous variables

(cid:26)

n

Ct, Yt, Πt,

Ci,t, Li,t, Li,t, Ui,t, Vi,t, Ni,t, Yi,t, {Xij,t}J

j=1

Given paths for the exogenous processes,

(cid:8)

Zt, εm
t

{αi,t, ri,t, si,t, χi,t, W s
i,t

1. Firms maximize profits,

2. Households maximize lifetime utility,

3. Markets clear

JX

(cid:27)∞

o

J

i=1

t=0
(cid:9)∞
t=0

}J

i=1

, such that:

Yi,t = Ci,t +

Xji,t,

Ni,t + ri,tVi,t = Li,t,

Yt = Ct

(30)

j=1

4. The central bank sets the nominal policy rate subject to the Taylor rule in Equation

(26).

11I also impose that

P

changing aggregate demand.

J
i=1 αi,t = 1 for all t so that αi,t shifts relative preferences across sectors without

16

3 Results: Inflation and labor market dynamics

In this section, I solve the model nonlinearly under perfect foresight and highlight three key
predictions. First, accounting for labor market frictions alters the propagation of sectoral
shocks. Second, the model features a nonlinear Phillips curve at the sectoral and aggregate
levels that steepens at high levels of labor market tightness. Third, as a result, how important
each sector is for aggregate inflation changes endogenously with labor market conditions. As
a result, monetary policy is less effective at stimulating output when some labor markets
are tight.
I then show how these three features together imply that a shift in relative
demand from one sector to another can generate aggregate inflation, a steepening of the
aggregate Phillips curve, and a decline in aggregate matching eﬀiciency—an outward shift
in the Beveridge curve.

3.1 Solution Method and Calibration

I solve the model nonlinearly by simulating impulse responses to a range of shocks under per-
fect foresight using the Newton-Raphson method. That is, I write the equilibrium conditions
as a system of nonlinear equations.

F (X) = 0,

(31)

where X stacks all endogenous variables over the simulation horizon T .12 I then solve for X
such that F (X) = 0 holds given the path of exogenous variables, using the steady-state as
the initial guess and the Jacobian to update the guess until convergence.

I solve the model nonlinearly for two reasons. First, search-and-matching introduces
potentially interesting nonlinear dynamics for large changes in tightness that I want to be
able to capture. Second, sector-level shocks can lead to large deviations from steady-state at
the sector level, even when aggregate variables remain relatively close to steady-state. The
large deviations from steady-state at the sector level necessitate a treatment that accounts
for the possible nonlinearities for large shocks. Labor market tightness for instance routinely
rises and falls by large amounts over the business cycle, and rose to over 100% above the
pre-Covid average in the aftermath of the pandemic.

12I use T = 200 and find that the results are robust to increasing the horizon T .

17

Figure 1: Fit of model implied Beveridge curve to JOLTS data (2000–2025)

Note: Orange dots plot monthly data on model consistent measure of aggregate unattached
workers (U +H in CPS/JOLTS data) and aggregate vacancies (U +H in JOLTS data). These
are the model constent estimates since JOLTs and CPS measures of V and U are vacancies
and unemployed worker who still have not matched at the end of the sample period, whereas
H is total hires over the entire sample period. The start of sample period U and V in
model terms are therefore U + H and V + H. The purple line plots an example Beveridge
curve implied by the model specification at the average level of hires. In reality, the model
beveridge curve is constantly shifting in response to shocks, and is not in fact stable at the
drawn curve. The purpose is to demonstrate that the shape of the curve (steeper when labor
markets are tight) conforms with the data.

I demonstrate the model’s mechanisms in a simulated symmetric version of the model
where each sector is calibrated to look like an average sector in the U.S. economy. I calibrate
the model so that each sector has an intermediate-share of 50%, which roughly matches the
intermediates share in the Bureau of Economic Analysis (BEA) input–output (I–O) tables.
I assume the share of their own output in intermediates is 70% and the share of the other
sector’s output in their intermediates is 30%, to reflect the large diagonal elements in the
I–O tables.

I set the other parameters of the model using existing estimates from the literature.
I set the elasticity of substitution between intermediates and labor (ϵy) to 0.2 to match
the estimated elasticity of substitution between intermediate inputs in Atalay (2017).
I
purposefully choose an estimate on the low end of the range of estimates in the litarature
because this parameter captures the elasticity of substitution between intermediates and
I want the model to
labor in the short run (a period in the model is just one month).

18

Unemployment(U)0.050.10.150.2Vacancies(V)0.040.060.080.10.120.14AggregateBeveridgeCurve,H=0.039ExampleModelBeveridgeCurveDatacapture that labor scarcity in a sector can lead to supply constraints because firms find it
diﬀicult to substitute between inputs in the short run, a feature highlighted as potentially
important in the post-pandemic recovery by Lorenzoni and Werning (2024). I set the sectoral
elasticity of substitution between intermediate suppliers (ϵ) to 4.33, which implies a markup
of about 30% (R. Hall, 2018).
I set ϕπ, ϕy, and ρi to 1.39, 1.01, and 0.82, respectively,
based on estimates for the Greenspan–Bernanke era in C. Carvalho et al. (2021). I set ψp
by aggregating the estimates for the frequency of price changes at the sector level in Pasten
et al., 2020 using sales shares as weights. I set φ = 3.57, consistent with Chetty (2012). I set
ψL = 31.24, roughly half the value of a job to the household, to capture the fact that workers
face high costs of switching across occupations and sectors (Artuç et al., 2010; Caliendo et
al., 2019; Cardoza et al., 2022; Humlum, 2021).

Table 1: Two-symmetric-sectors Calibration

Parameter Value Source

ϵd
ϵy
ϵ
ϕπ
ϕy
ρi
κ
ρw
σ
φ
ηi
Ωd,i
Ωx,i
Ωii
Ωij
ψp,i
ψL,i

Implied markup of 30%, R. Hall (2018)

0.60 Atalay (2017)
0.20 Elasticity of substitution between intermediate inputs, Atalay (2017)
4.33
1.39 Estimates for Greenspan–Bernanke Period, C. Carvalho et al. (2021)
1.01 Estimates for Greenspan–Bernanke Period, C. Carvalho et al. (2021)
0.82 Estimates for Greenspan–Bernanke Period, C. Carvalho et al. (2021)

Equal surplus sharing between households and firms

1
0.8 —
2 —

3.57 Chetty (2012)
0.93 Estimate using JOLTS (2000-2025)
0.50 Nominal Consumption Share (Symmetric)
0.50 Roughly Intermediates Share in BEA I–O tables
0.70 —
0.30 —
60.11 Pasten et al. (2020)
31.24 Half value of a job to household.

Note: Calibrated parameter values in two-sector-symmetric economy. Based on estimates or standard values
in the literature.

I assume a matching function of the form mi,t(Ui,t, Vi,t) =

ηi , which sat-
isfies the properties outlined in Section 2. Crucially, it ensures that Hi,t ≤ min{Ui,t, Vi,t}
while remaining differentiable throughout.
I estimate ηi by nonlinear least squares using
aggregate data on vacancies, unemployment, and hires from the Job Openings and Labor
Turnover Survey (JOLTS) from 2000–2025. Figure 1 demonstrates that the model-implied
Beveridge curve at the average hiring rate is a good fit to the JOLTS data. In particular, the

U

(cid:0)

(cid:1)− 1

−ηi
i,t + V

−ηi
i,t

19

functional form I assume can capture both the flatter portion of the Beveridge curve when
unemployment is high, and the steepening at low levels of unemployment. I summarize the
full parametrization in Table 1.

3.2 The Propagation of Shocks Through the Labor Market

I begin by showing how the presence of labor market frictions alters the propagation of
sectoral shocks via a new labor channel. Higher labor demand in one sector leads to a
reduction in labor supply elsewhere, raising hiring costs across the network. As a result,
positive demand shocks in one sector can effectively create negative supply shocks in other
sectors. This is the converse of, but closely related to, the finding in Guerrieri et al. (2022)
that negative supply shocks in one sector lead to negative demand shocks in other sectors.
Consider, for instance, a permanent 1% increase in the household consumption preference
i,t. The shock raises αi,t, and

for goods in sector i, captured by a permanent increase in εα
therefore, as Figure 2a shows, household consumption of sector i’s output, Ci,t.

Ci,t = αi,t

(cid:19)−ϵd

(cid:18)

Pi,t
Pt

Ct,

Cj,t = αj,t

(cid:19)−ϵd

(cid:18)

Pj,t
Pt

Ct

P

J

j=1 αj,t = 1, this implies a decline in αj,t for some j ̸= i, and therefore a decline in
Because
household consumption of sector j’s output.13 For simplicity, I assume in the figures below
that the increase in αi,t is offset entirely by an equal decline in αj,t in just one other sector j.
The increase in demand for sector i output leads to an increase in inflation in sector i
(see Figure 2b). As is well understood, in any production network economy this increase in
prices spills over to other sectors through input–output linkages in production. The increase
in inflation in sector i affects the price of intermediate inputs, and therefore marginal costs,
in sector j. All else equal, this results in higher sector j inflation than absent the network,
raising aggregate inflation.

As firms in sector i increase hiring to meet the increase in consumption demand for their
good, both wages and the job-finding rate in sector i rise (see Figure 3b). The job-finding rate
increases on impact because the increase in vacancy postings by firms raises labor market
tightness, θi,t, and therefore qi,t. Because the households’ optimal search strategy implies
that the amount of labor searching in sector i rises with wages and the job-finding rate, this
in turn, leads to an increase in the labor supply to sector i (See Figure 2c), and, similarly, a
decline in labor supply in sector j.

13As a result, this shock is pure demand-reallocation shock that does not effect the level of aggregate

demand.

20

Figure 2: Impulse Responses to an increase in αi,t

(a) Consumption & shock in i

(b) Inflation in i

(c) Labor Supply in i

(d) Tightness in i and j

Note: Impulse responses to an increase in relative household preferences for consumption in sector
i. Top left panel: Impulse response of real consumption in sector i Ci,t (in purple) to the preference
shock in sector i αi,t (in orange). Top right panel: Impulse response of inflation in sector i Πi,t.
Bottom left panel: Impulse response of labor supply in sector i Li,t. Bottom right panel: Impulse
response of tightness in sector i θi,t (in purple) and in sector j θj,t (in orange). All variables in
percentage deviations from steady-state.

The spillovers in labor supply from sector j to sector i increase tightness, and therefore
hiring costs and marginal costs in sector j, leading to a new spillover from sector-specific
shocks. Indeed, as Figure 2d demonstrates, the labor supply spillovers can be strong enough
to eventually fully offset the increase in labor demand in sector i, leading to a decline in
tightness in sector i and an increase in tightness in sector j, once the initial spike in labor
demand subsides in sector i.

21

Period051015202530%dev.froms.s.00.20.40.60.811.2ConsumptionPreferenceShockPeriod051015202530%dev.froms.s.0123456Period051015202530%dev.froms.s.00.10.20.30.40.50.6Period051015202530%dev.froms.s.-20-100102030SectoriSectorjFigure 3: Impulse Responses to an increase in αi,t

(a) Employment & vacancies in i

(b) Wages & job-finding rate in i

Note: Impulse responses to an increase in relative household preferences for consumption in sector
i. Left panel: Impulse response of employment in sector i Ni,t (in purple, left y-axis) and vacancies
sector i Vi,t (in orange, right y-axis). Right panel: Impulse response of real wages in sector i Wi,t
Pt
(in purple, left y-axis) and job-finding rate in sector i Fi,t (in orange, right y-axis). All variables in
percentage deviations from steady-state.

3.3 A nonlinear Sector-Specific Phillips Curve

The presence of labor market frictions leads to a nonlinear Phillips curve that steepens at
high labor market tightness. This finding is consistent with the evidence in Gitti (2024) and
Benigno and Eggertsson (2023). Benigno and Eggertsson (2023) develop a model to capture
this nonlinearity by assuming a kink in the wage-setting process at levels of tightness above 1.
I demonstrate that the nonlinearity arises from a standard search-and-matching framework
with generalized Nash bargained wages once I solve the model nonlinearly. In addition, in a
multi-sector economy, the presence of nonlinearities at the sector level leads to endogenous
changes in which sectors are the most important for aggregate inflation as local labor market
conditions change.

To build intuition, consider a simplified version of the model outlined above where wages
are rigid, ρw = 1, there are no household labor adjustment costs, ψL = 0, and firms and
households make static labor supply decisions. This last simplification can be rationalized
by assuming, as in Benigno and Eggertsson (2023), that all workers separate at the start of
each period, before a random fraction (1 − si,t) are reemployed in sector i without needing
to go through the matching process. As a result, households and firms no longer account
for the dynamic consequences of their labor supply and demand decisions. In this case, to

22

Period0102030ProductionWorkers00.10.20.30.4VacancyPostings05101520Employment(left)Vacancies(right)Period0102030Realwages00.511.5Job--ndingrate-20246810Wages(left)Job--nding(right)first-order aggregate inflation is given by

πagg
t = Γθθt + ΓπEtπt+1 + υt.

(32)

Where θt is a J × 1 vector of sectoral tightness, πt+1 is a J × 1 vector of sectoral inflation,
and υt is an endogenous cost push shock as in Rubbo (2023). I derive this first-order ap-
proximation of the aggregate New Keynesian Phillips curve by combining Equations (12)
and (16). It holds for any constant returns to scale production and matching functions, and
does not depend on the specific functional forms assumed above. 14 Γθ and Γπ are 1 × J
coeﬀicient vectors that capture how important each sector is for aggregate inflation. This
expression is equivalent to the network Phillips curve derived in Rubbo (2023), adjusted for
the presence of labor market frictions.

The coeﬀicient matrices on tightness and expected future inflation are

Γθ = Ω′
Γπ = βΩ′

d [I − λ (Ωx − I) (I − 1Ω′
d)]
d [I − λ (Ωx − I) (I − 1Ω′

−1 ΩnΓQη
−1
d)]

where Ωd is a J × 1 vector of steady-state nominal consumption shares, Ωx is the J × J
input–output matrix capturing the steady-state expenditure shares of each sector i on each
intermediate input j, Ωn is a J × J diagonal matrix with the steady-state labor share in each
sector on the diagonal. The two new terms relative to a standard production network setup
are η, a J ×J diagonal matrix with the negative of the elasticity of the vacancy-filling rate to
changes in tightness along the diagonal, and ΓQ, a J × J matrix capturing the pass-through
from changes in the vacancy-filling rate to marginal costs in each sector.

The presence of labor market frictions therefore alters the propagation of sectoral shocks
to aggregate inflation in two distinct ways. First, as demonstrated in the previous section, an
increase in labor demand in one sector leads to a reduction in labor supply in other sectors,
potentially leading to a cascade of labor market tightness throughout the network. As a
result, sectors that lead to larger spillovers in tightness to other sectors, and therefore to
larger changes in the entire vector of sectoral tightness, θt, are more important for aggregate
inflation. Second, how changes in tightness pass through to prices depends on the details of

14See the appendix for a detailed derivation. The expression for the endogenous cost push shock is

υt = Ω′
+ Ω′
+ Ω′

d [I − λ (Ωx − I) (I − 1Ω′
d [I − λ (Ωx − I) (I − 1Ω′
d [I − λ (Ωx − I) (I − 1Ω′

d)]
d)]
d)]

−1 [I + λ (Ωx − I) (I − 1Ω′
−1 [ΩnΓQrt − at]
−1 [I − λ (Ωx − I) 1Ω′

d] αt

d)] pt−1

23

the sector-specific matching process. Labor market conditions in sectors where the vacancy-
filling rate changes more in response to changes in tightness, where hiring costs respond more
to changes in the vacancy-filling rate, and where marginal costs respond more to changes in
labor costs, matter more for aggregate inflation.

To first order, shifting demand into sectors where labor markets are more rigid—where
firms have a harder time adjusting their labor input either because increasing vacancy post-
ings leads to fewer additional hires or because the effective hiring costs rise more quickly
as tightness increases—leads to an endogenous cost-push shock, triggering aggregate infla-
tion. The intution is similar to Rubbo (2024), who shows that shifting demand into sectors
with less elastically supplied inputs leads to higher aggregate inflation. Intuitively, sectors
with less elastically supplied inputs have a harder time increasing output in response to a
positive demand shock, leading to a smaller increase in quantities and a larger increase in
prices when these sectors experience a surge in demand. In Rubbo (2024), all sectors face
identically elastic labor supply curves, but vary in terms of their capital and intermediate
input usage, leading to differences in the elasticity of input supply curves across sectors.
Here, I show that a similar mechanism arises when sectors face frictional labor markets: any
variation in how rigid labor markets are across sectors leads to differences in the elasticity of
input supply curves, and therefore alters how important each sector is for aggregate inflation.
In the appendix, I demonstrate using the BLS’s Job Openings and Labor Turnover Survey
(JOLTS) data that there is substantial heterogeneity in average vacancy-filling rates across
sectors, and therefore in steady-state hiring costs.15

In the nonlinear full model, the elasticity of the vacancy-filling rate and the pass-through
to hirings costs changes endogenously with local labor market conditions. Figure 4a demon-
strates that as tightness in sector i rises, the vacancy-filling rate falls. As a result, firms
in sector i need to post more vacancies to achieve the same level of hiring:
it gets harder
and harder for firms in that sector to hire additional workers. This, in turn, leads to both
a direct increase in marginal costs (see Figure 4b), operating both through hiring costs and
wages, and reduces the elasticity of the labor supply curve faced by firms in sector i. Both
the hiring costs and the elasticity of labor supply change nonlinearly as the constraint that
Hi,t ≤ mi,t(Ui,t, Vi,t) approaches, leading to endogenous changes in which sectors are the
most important for aggregate inflation.

15The average monthly vacancy-filling rate in the finance and insurance industry, for instance, from 2000
to 2025 is about 0.37, while the average vacancy-filling rate in construction over the same horizon is about
0.68. These large differences in average vacancy-filling rates likely reflect differences in how search works
across sectors, depending on the skill level of employees and the time each application takes to process.

24

Figure 4: How Hiring costs and Marginal Costs vary with Tightness

(a) Relationship between θi,t and qi,t

(b) Relationship between θi,t and M Ci,t

Note: Right panel: Relationship between tightness and vacancy-filling rate in sector i, gener-
ated by varying size of demand shock in sector i. Left panel: Relationship between tightness
and marginal costs in sector i, generated by varying size of demand shock in sector i. All
variables in percentage deviations from steady-state.

As tightness rises, firms become effectively more supply constrained. As a result, they
change output by less and inflation by more in response to additional demand or supply
shocks. Indeed, as Figure 5 demonstrates, this mechanism results in a nonlinear Phillips
curve at the sector level. The solid purple line in Figure 5 plots the relationship between
inflation and unemployment in response to sector-specific demand shocks. The inflation-
unemployment Phillips curve steepens as unemployment falls precisely because firms become
more supply constrained at high levels of tightness, forcing more of the effects of positive
demand shocks into prices rather than quantities.

Relatedly, because output becomes more constrained at high levels of tightness, monetary
policy has larger effects on inflation and smaller effects on output when tightness is high. To
demonstrate this effect, Figure 6a plots ∂Π
∂ipol (θ), the effect of a change in the nominal policy
rate on inflation, by plotting the response of inflation to a 1 percentage point policy rate
cut, relative to the effect of a rate cut when tightness is initially at steady-state. Conversely,
Figure 6b plots ∂Y
∂i (θ), the effect of a change in the nominal policy rate on output, by plotting
the response of output to a 1 percentage point policy rate cut, relative to the effect of a rate
cut when tightness is initially at steady-state. When tightness in sector i is 125% above
steady-state, for instance, monetary policy has a nearly twice as large effect on inflation in
that sector.

25

Tightness(3i;t)-50050100Vacancy-FillingRate(qi;t)-40-30-20-100102030Tightness(3i;t)-50050100MarginalCost(MCi;t=Pt)-10-505101520Figure 5: Nonlinear sector-specific Phillips curve in unemployment

Note: Phillips curve relationship between inflation and unemployment in sector i, generated
by varying size of demand shock in sector i. The curve steepens at high levels of tightness,
and therefore low levels of unemployment, because firms become more supply constrained
as tightness rises.

The nonlinear Phillips curve is consistent with the empirical findings in Gitti (2024) that
the Phillips curve is steeper in U.S. in regions with tighter labor markets. In addition, as
argued in Benigno and Eggertsson (2023), the presence of a nonlinear Phillips curve can help
explain why central banks were surprised by the surge in inflation in 2021. Assuming the
Phillips curve is linear curve underestimates the inflationary pressure in tight labor markets,
and leading to larger than expected inflationary pressure from aggregate demand stimulus.
However, my findings are distinct in two important ways. First, unlike Benigno and Eg-
gertsson (2023), where the Phillips curve steepens when tightness crosses a certain threshold
because of an assumed kink in the wage setting process, I highlight an alternative channel
generates a non-linear Phillips curve. The Phillips curve steepens in the standard search-
and-matching framework once I solve the model nonlinearly because search-and-mathing
frictions essentially introduce capcity constraints in the labor market. Firms in tight labor
markets struggle to hire, and therefore to increase output as much and as quickly as they
would like, and instead push additional changes in demand or supply into prices rather than

26

UnemploymentUi;t-30-20-100102030In.ation(&i;t)-2-101234quantities. This suggests that the findings in Benigno and Eggertsson (2023) do not rely
on the specific specification for the wage setting and search process, and instead arise more
generally in models featuring frictional labor markets.

Figure 6: Effects of monetary policy on output and inflation in sector i

(a) Inflation in sector i

(b) Output in sector i

Note: Right panel: Effect of -1% monetary policy shock on inflation in sector i at different
levels of tightness in sector i, relative to effect of monetary policy at steady state. Left panel:
Effect of -1% monetary policy shock on Ci at different levels of tightness in sector i, relative
to effect of monetary policy at steady state.

Second, by solving the model nonlinearly instead of assuming a kink, I show that search-
and-matching frictions lead to a continuous steepening of the Phillips curve: the curve gets
progressively steeper as tightness rises. As a result, assuming a linear Phillips curve becomes
a worse approximation as tightness rises substantially above steady-state, as it did in the
U.S. in 2021 and 2022. In these circumstances, the monetary authority is especially prone
to allowing inflation to surge more than expected if it extrapolates from the relatively flat
Phillips curve in normal times.

3.4 Aggregate Effects of Sectoral Demand Shocks

I now demonstrate how a production network, capturing sector-specific shocks, and labor
market frictions alter the aggregate effects of demand shocks.
In particular, I show that
the Phillips curve steepens and monetary policy becomes less effective at stimulating output
when just a few sectors experience large increases in tightness.

As demonstrated above, a shock to relative demand in sector i, αi,t, increases demand in
sector i and reduces demand in sector j. As the size of the relative demand shock increases,

27

Tightness(3i;t)-40-20020406080100120Relativetoe,ectats.s.-20020406080100Tightness(3i;t)-40-20020406080100120Relativetoe,ectats.s.-30-25-20-15-10-50510so does the increase in labor demand and tightness in sector i. For instance, in the baseline
calibration, an 8% increase in αi,t triggers about a 180% increase in tightness.

Figure 7: Aggregate effects of sectoral demand shocks

(a) Tightness in i

(b) Aggregate inflation

(c) Matching eﬀiciency

Note: All panels plot response to αi,t shocks of various sizes. Right panel: Effect of shifting
consumer preferences on tightness in sector i. Middle panel: Effect of shifting consumer
preferences on aggregate inflation. Left panel: Effect of shifting consumer preferences on
measured aggregate matching eﬀiciency Ht√
.

UtVt

As tightness rises in sector i, it gets progressively more diﬀicult to hire, leading to both an
increase in inflation and a decrease in aggregate matching eﬀiciency, a measure of how many
hires are generated for a given number of unemployed workers and vacancies. Aggregate
inflation rises because the rise in tightness in sector i causes firms in that sector to become
relatively more supply constrained. Firms in sector j, on the other hand, experience a decline
in tightness. As a result, prices rise by more in sector i than they fall in sector j. In addition,
the rise in αi,t gives sector i a larger weight in the aggregate price index. These two factors
combine to produce a significant rise in aggregate inflation:
in the illustrative symmetric
calibration an 8% increase αi,t, in line with the shift from goods to services following the
Covid pandemic, generates about a 7% increase in aggregate inflation.

Figure 7c demonstrates that a shift in demand from one sector to the other also results in
. Aggregate

a decline in the aggregate matching eﬀiciency, which I define as ϕagg

t =

H agg
t√
U agg
t

V agg
t

matching eﬀiciency falls because labor demand is concentrated precisely in the sector where
it is hardest to hire. It gets harder to hire in sector i for two reasons. First, as in Şahin et al.
(2014), sectoral shocks can generate mismatch between where vacancy postings are and where
unattached workers search. This mismatch is exacerbated when labor adjustment costs are
larger. Second, the number of hires generated per additional vacancy declines endogenously
as the constraint Hi,t ≤ mi,t(Ui,t, Vi,t) approaches. That is, the natural constraint that firms
cannot create additional workers out of thin air by posting a larger number of vacancies,

28

Sizeofshock-505%dev.froms.s.-1000100200Sizeofshock-505%devfroms.s.0246Sizeofshock-505%devfroms.s.-20-100naturally implies that a large increase in vacancies in just one sectoral labor market leads
to a decline in the number of hires relative to what one would expect given the number
of aggregate vacancies and unemployed workers. An 8% increase αi,t generates about a
25% decrease in measured aggregate matching eﬀiciency, despite no change in underlying
matching eﬀiciency at the sector level.

Figure 8: Aggregate effects of monetary policy on output and inflation

(a) Aggregate inflation

(b) Aggregate output

Note: Right panel: Effect of -1% monetary policy shock on aggregate inflation relative to
effect at steady state. Purple line is effect in model with labor market frictions, orange line
is effect in model without labor market frictions. Left panel: Effect of -1% monetary policy
shock on aggregate output relative to effect at steady state. Purple line is effect in model
with labor market frictions, orange line is effect in model without labor market frictions.

As figure 8 demonstrates, in the presence of labor market frictions, monetary policy has
smaller effects on output and larger effects on inflation when just a few sectors experience
larger increases in tightness. The solid purple line in the top panel plots the effect of a 1
percentage point cut in the nominal policy rate on aggregate inflation, relative to the effect
of a rate cut when αi,t is at steady-state. Conversely, the solid black line in the bottom
panel plots the effect of a 1 percentage point cut in the nominal policy rate on aggregate
real output, relative to the effect of a rate cut when αi,t is at steady-state. Following an 8%
increase in αi,t, monetary policy has about 160% larger effects on aggregate inflation than it
does at the steady state.

As tightness rises in sector i, and firms in that sector become more supply constrained,
they become progressively less able to increase output in response to an additional positive
aggregate demand shock, and therefore must raise prices by more instead. As a result, it is

29

SizeofShockto,i;t-8-6-4-202468Relativetoe,ectats.s.020406080100120140160180WithLaborMarketFrictionsWithoutLaborMarketFrictionsSizeofShockto,i;t-505Relativetoe,ectats.s.-40-30-20-100WithLaborMarketFrictionsWithoutLaborMarketFrictionsenough for some sectors to be supply constrained for the aggregate Phillips curve to steepen
and for monetary policy to become less effective at stimulating output. This suggests that
central banks must account for sectoral labor market conditions when predicting the effects
of their monetary policy actions. When there are exceptional labor market conditions in just
a few sectors in the economy, acting as if the aggregate Phillips curve remains flat may lead
the central bank to underestimate the inflationary consequences of its actions. Conversely, in
an uneven downturn or recovery, monetary policy may be a less effective tool for stimulating
aggregate output since constrained sectors respond less to additional demand stimulus.

The last prediction is an important insight for when monetary policy is and is not likely
to be an effective stabilization tool. For instance, consider two shocks that lead to an identi-
cal increase in aggregate demand. Suppose the first shock is a true aggregate demand shock
affecting all sectors equally, while the second raises demand in some sectors by significantly
more than others. The model’s predictions for the effects of sectoral demand shocks suggest
that monetary policy will have a larger effect on output in the case of the pure aggregate
demand shock, where no sector experiences a spike in tightness and becomes supply con-
strained, than in the case of an uneven sectoral demand shock. Similarly, if sectoral shocks
lead an overall decline in demand, triggered by large declines in demand in some sectors
and an increase in demand in others, monetary policy may be unable to stimulate aggregate
output, and stimulative monetary policy will instead have larger inflationary effects.

4 Quantitative Application to Post-Pandemic Inflation

and Labor Market Dynamics

I now turn to a quantitative application of the model to the recent inflation and labor market
dynamics in the United States following the COVID-19 pandemic. I begin by documenting
sectoral differences in inflation and labor market dynamics during the recovery.
I then
demonstrate that a shock to relative demand for goods over services can partially account
for both inflation and matching eﬀiciency in an estimated two-sector goods-services version
of the model.

4.1 Inflation and Labor Market Dynamics in Goods and Services

From 2021 to 2024, the United States experienced two phenomena unprecedented since at
least the 1980s: (1) inflation surged to levels not seen since the Volcker disinflation, and (2)
labor market tightness reached historic highs, with job openings per unemployed worker at
levels unseen since the 1960s. I begin by documenting sectoral differences in these dynamics

30

during the pandemic recovery, focusing on the responses in the goods and services sectors.
Inflation, real consumption, labor market tightness, and new-hire wages rose more in goods
than in services, suggesting a link between sectoral inflation and labor market dynamics.

Figure 9: Labor Market and Inflation Dynamics in Goods and Services

(a) PCE Inflation.

(b) Real PCE Shares

(c) Tightness

(d) New Hire Wages

(e) Matching Eﬀiciency

Note: The gray shaded area indicates the period from the start of the inflation surge in early
2021 to the first Federal Reserve rate hike in March 2022. Top left: PCE inflation in goods
(purple) and services (orange). Top middle: Real PCE consumption share of goods (purple)
and services (orange), normalized to 100 in January 2020. Top right: Labor market tightness
in goods (purple) and services (orange), normalized to 100 in January 2020. Bottom left:
Real new hire wages in goods (purple) and services (orange), detrended from 2000 to 2025,
normalized to 100 in January 2021 to demonstrate change in real wages over inflation period.
Bottom right: Aggregate matching eﬀiciency, measured as the ratio of hires to the square
root of the product of unemployed workers and vacancies, normalized to 100 in January
2020.

In the figures above, I use the PCE goods and services price and quantity indices from
the Bureau of Economic Analysis. I use vacancy data from the Bureau of Labor Statistics’
Job Openings and Labor Turnover Survey (JOLTS) roughly at the 2-digit NAICS level,
and aggregate to the two-sector goods-services breakdown available from the PCE. I use
unemployment data from the Current Population Survey (CPS).

As Figure 9 shows, the surge in U.S. inflation, which began in early 2021. Year-over-
year inflation, which peaked at 8.9 percent in June 2022, was uneven across the economy.

31

2019202020212022202320242025Date−2.50.02.55.07.510.0PCEInﬂation(YoY%)GoodsServices2019202020212022202320242025Date95100105110115Jan2020=100GoodsServices2019202020212022202320242025Date50100150200250300Tightness(Jan2020=100)GoodsServices2019202020212022202320242025Date98100102104106DetrendedRealWagesGoodsServices2019202020212022202320242025Date8090100MatchingEfﬁciencyPCE inflation in goods jumped from about 0 to over 10 percent by mid-2022, while services
inflation rose more slowly, peaked around 6 percent, and stayed persistently elevated near
4 percent into 2024. This pattern suggests that a mechanism where inflation in one sector,
triggered by either sectoral supply or demand shocks, gradually propagates to other sectors
through the production network as in Minton and Wheaton (2023).

This initial surge followed a sharp rise in demand for goods, as consumers shifted toward
lower-contact spending. The goods share of real PCE rose from 32 to 37 percent early
in the pandemic—a 15 percent jump—and has remained elevated since. Alongside this
shift, supply chain disruptions and the war in Ukraine created supply shocks that hit goods
disproportionately.16 Recent work emphasizes the importance of these sector-specific demand
and supply disturbances for the inflation surge (Amiti et al., 2023; Comin et al., 2023; Di
Giovanni et al., 2023; di Giovanni et al., 2023; Ferrante et al., 2023; Guerrieri et al., 2022;
Lorenzoni & Werning, 2024; Rubbo, 2024).

Figure 9 shows that inflation rose alongside labor market tightness, measured as the
I calculate searchers using CPS
sector-specific ratio of job vacancies to total searchers.
microdata on unemployment-to-employment (U E), nonparticipation-to-employment (N E),
and employment-to-employment (EE) transitions, adjusted following Fujita et al. (2024).
Assuming random matching and sector-specific job search, total searchers in sector i are:

T Si,t = Ui,t +

ρN E
i,t
ρU E
i,t

Ui,t
Ut

Et +

ρN E
i,t
ρU E
i,t

Ui,t
Ut

Nt

(33)

where Ui,t is the number of unemployed workers most recently employed in sector i, Ut is
the total number of unemployed workers, Et is the total number of employed workers, Nt is
the total number of people not currently in the labor force, ρN E
is the N E rate into

i,t =

H N
i,t
Nt

H E
i,t
Ut

is the U E

sector i (H N

i,t is the number of hires in sector i from nonparticipation), ρU E

i,t =

H E
i,t
Et

i,t =

rate into sector i, and ρEE

is the EE rate into sector i.
I use this broader measure rather than the conventional unemployment-based one for
three reasons. First, recent work highlights the role of EE transitions in post-pandemic
labor markets (Autor et al., 2023; Bagga et al., 2025; Barlevy et al., 2023; Faccini & Melosi,
2025; Moscarini & Postel-Vinay, 2023). Second, Barnichon and Shapiro (2024) show that
tightness based on total searchers forecasts inflation better than unemployment alone. Third,
the measure aligns with the model introduced above, in which households allocate members
to search across labor markets, and some unattached members find a new job within the

16For example, the New York Fed’s Global Supply Chain Pressure Index peaked at 4.5 standard deviations
above average in December 2021, the highest on record. The ISM Supplier Deliveries Index reached a record
78.8 in May 2021.

32

same period, never showing up in end-of-period unemployment numbers.

The goods-sector labor market was about 80 percent tighter in early 2022 than in January
2020, when conditions were already historically tight.17 Services were about 40 percent
tighter. Goods-sector tightness closely tracks goods inflation, underscoring the potential
importance of accounting for labor market dynamics, even in a multi-sector setting.

Similarly, Figure 9e shows that matching eﬀiciency—the number of hires per vacancy and
searcher—declined sharply, and by a comparable amount to the Great Recession. Figure 9e
plots the residuals from a regression based on a Cobb-Douglas matching function,

Hi,t = ϕi,tT Sη

i,tV 1−η

i,t

(34)

where η is the elasticity of the matching process to total searchers and ϕi,t is the matching
eﬀiciency.

I run the following regression separately for the aggregate data:

(cid:19)

(cid:18)

Ht
T St

log

= log ϕ + (1 − η) log θt + ϵt

(35)

where log ϕ + ϵt = log ϕt is the log matching eﬀiciency estimated in the data. Measured
aggregate matching eﬀiciency declines substantially following the pandemic. My model links
this fall in matching eﬀiciency, which shifts the Beveridge curve, to higher firm marginal
costs and inflation.

Despite the historic rise in labor market tightness and the substantial decline in matching
eﬀiciency, a common objection to labor-based explanations is that real wages fell as prices
Indeed, aggregate real wages declined in both goods and services. But Figure 9d
rose.
shows that real wages for new hires—relevant for marginal costs—rose sharply just as goods
inflation accelerated. Because services dominate employment, the aggregate masks these
sectoral patterns. As the model demonstrates, even without wage growth, rising tightness
and falling matching eﬀiciency raise costs by increasing firms’ hiring costs.

Rising tightness also reduces the effective elasticity of labor supply to firms, and especially
so in the sector experiencing a larger rise in tightness. Cross-sectoral shifts in tightness driven
by changes in relative demand can therefore lead to an endogenous decrease in the input
elasticity of the goods sector, exerting upward pressure on prices even absent a wage response.
Taken together, these facts point to both sectoral heterogeneity and labor market frictions
as central to the pandemic recovery, consistent with the model outlined above.

17I show in Appendix C that the same pattern holds when using the conventional V /U tightness measure.
Tightness moves more, though, when using V /U , and by this metric the goods-sector labor market was about
150 percent tighter in early 2022 than in January 2020.

33

4.2 Can Demand Shocks Explain the Joint Dynamics of Inflation

and the Labor Market?

As I demonstrate above, the model’s predictions are broadly consistent with the experience
during the post-pandemic recovery.
In this section, I use Bayesian methods to estimate
parameters in a linearized version of the model on data from 2000 to 2019, a period when
shocks were small relative to the COVID-19 period and where the linear approximation is
therefore more likely to be accurate. I then use these estimated parameters to calibrate the
fully nonlinear model to assess the impacts of different shocks during the post-pandemic
recovery. I focus on peak responses since my model lacks many of the mechanisms used in
the literature to smooth out responses over time.

The first three columns of Table 2 report the prior mean, variance, and distribution for
each of the estimated parameters. The last column reports the mode from 3 million draws
from the posterior generated with a standard Random-Walk Metropolis–Hastings algorithm.
The priors are taken from the literature and are listed in Table 1.

I then ask whether the model can quantitatively account for changes in inflation and labor
market conditions with shocks to consumer preferences for relative consumption and then
with a combination of preference shocks, aggregate demand shocks, aggregate labor supply
shocks, and a positive shock to the separation rate. I find that the shift in consumption
demand from services to goods alone can account for about 2.9 percentage points of inflation,
roughly 40% of the total increase in inflation. The shift in consumption demand can also
broadly account for movements in aggregate matching eﬀiciency and relative changes in
tightness in goods and services. It cannot, however, account for the magnitude of changes
in tightness and counterfactually predicts a fall in tightness in services.

Allowing for additional declines in aggregate labor supply and increases in separation
rates, consistent the post-Pandemic labor market data, and an increase in aggregate demand,
allows the model to match the magnitude of changes in tightness and inflation.18 The
accounting for the shift in consumption

18Afrouzi et al. (2024) provide one possible rationale for the observed increase in separation rates, positing
that workers are more likely to engage in on-the-job search and accept outside offers when inflation erodes
the value of their existing wages. Bagga et al. (2025) provide an alternative explanation based on shifting
preferences for job amenities during the pandemic.

34

Table 2: Prior Distributions and Estimated Posterior Parameter Values

Parameter Mean Variance Distribution

Mode

80% C.S.

ϕπ
ϕy
ρi
κ
σ
φ
ψL
ρw

1.39
1.01
0.82
0.50
2.00
3.57
36.17
0.80

0.30
0.10
0.15
0.10
0.50
0.50
10.00
0.15

Truncated Normal
Gamma
Beta
Gamma
Gamma
Gamma
Gamma
Beta

2.04
0.89
0.34
0.88
1.60
3.34
14.33
0.43

(1.80, 2.27)
(0.79, 1.02)
(0.20, 0.48)
(0.76, 1.00)
(0.43, 0.79)
(2.95, 3.82)
(11.31, 17.58)
(0.35, 0.51)

Note: The left panel shows prior distributions taken from the literature; see Table 1. Set
price rigidity parameter based on Pasten et al. (2020). Last column is the 80% credible set
from the posterior distribution.

I calibrate the remaining parameters of the model to match the real consumption share,
labor share, and input–output structure of the goods and services sectors, and to match
the sector level frequency of price changes in (Pasten et al., 2020). I report the network
parameter values in Table 3. I set the values of ϵd, ϵy, and ϵ to those reported in Table 1.

Table 3: Sectoral Shares and Elasticity Parameters

Category

Consumption Share
Intermediates Share (Total)
Intermediates Share (Goods)
Intermediates Share (Services)

Goods

0.3338
0.6210
0.4171
0.2039

Services

0.6662
0.3921
0.0684
0.3237

Notes: Goods–services sectoral consumption shares, and intermediate shares, calculated from PCE and
BEA input–output tables.

I start by introducing a shock to the relative demand for goods over services, αG,t,
calibrated to match the peak increase the real consumption share in goods post-pandemic
of 13.9 percent. I set the persistence of this shock to 0.98 to roughly match the path of the
real consumption share in goods from 2021 to the end of 2024. I then solve for the responses
to this shock under perfect foresight as described in the previous section.

35

Figure 10: Shock to Consumption Preferences Only

(a) Aggregate Inflation

(b) Labor Market Tightness in Goods

(c) Aggregate Matching Eﬀiciency

(d) Labor Market Tightness in Services

Note: All panels plot model predicted peak responses to shock to relative consumption
preferences for goods (in purple) and peak responses in the data (in orange). Red line is model
response at the posterior mode. Black dashed lines are at steady state levels. The purple
shaded areas are 68- and 95-percent credibility sets. Top left: Aggregate inflation. Orange
is peak CPI inflation of 8.9 percent in June 2022. Red line is model response at posterior
mode. Dashed black line is at steady state inflation rate of 2 percent. Top right: Increase
in tightness in goods sector. Bottom left: Peak decline in aggregate measured matching
V × Y ) in model and data (excluding early 2020 shutdown). Bottom right:
eﬀiciency (H/
Increase in tightness in services.

√

36

ModelDataPercentagePoints02468Mode95%c.s.68%c.s.ModelDataV=U00.511.5Mode95%c.s.68%c.s.ModelData%dev.froms.s.-10-50Mode95%c.s.68%c.s.ModelDataV=U00.511.5Mode95%c.s.68%c.s.The demand reallocation alone can account for about 2.9 percentage points of inflation,
or roughly between 40 percent of the total increase in inflation. Accounting for parameter
uncertainty, the shift in consumption demand can account for between 10 and 75 percent
of the inflation surge. While the reallocation shock alone cannot explain the entire inflation
surge, even at the upper end of estimates, this is hardly surprising given the many additional
disturbances that characterized the post-COVID period.

The demand reallocation can also broadly explain the relative patterns in tightness in
the goods and service sectors and the decline in aggregate matching eﬀiciency. The model
predicts that, as a result of increased relative demand in the goods sector, tightness in the
goods sector rises by between 45 and 57 percent more than in services (relative to steady state
tightness in both sectors). This accounts for a between 67 and 110 percent of the difference
in relative tightness observed in the data. In addition, the shift in demand from services to
goods explains between 70 and 93 percent of the decline in measured aggregate matching
eﬀiciency because of both missmatch between the sectors where vacancies are posted and
where unattached workers search, as well as the endogenous decline in hires per vacancy the
goods sector.

The shift in demand alone, however, cannot account for changes in the level of tightness
in the two sectors or in the aggregate. While tightness does rise substantially in the goods
sector, and the shock can explain approximatly 1/3 of the increase in tightness there, because
demand declines in services, tightness there only rises slightly due to the endogenous response
of households to search more in the goods sector. Since most employment is in services, the
shock explains only a small fraction of total aggregate tightness.

To give the model a better chance of matching the aggregate tightness I introduce two
additional shocks: (1) a persistent negative aggregate labor supply shock, calibrated to
match the roughly 3 percent decline in aggregate labor force participation observed during
the pandemic recovery, and (2) a persistent positive aggregate demand shock, calibrated to
match the remaining increase in inflation.

Together these shocks generate an increase in tightness in both the goods and service
sectors, that explains between 10 and 53 percent of the increase in tightness in services
and between 68 and 84 percent of the increase in tightness in the goods. Furthermore, the
model predicts that absent the shift in demand from services to goods, the aggregate demand
stimulus needed to match the rise in inflation would be about 50 percent higher.

37

Figure 11: Combination of Shocks

(a) Aggregate Inflation

(b) Labor Market Tightness in Goods

(c) Aggregate Matching Eﬀiciency

(d) Labor Market Tightness in Services

Note: All panels plot model predicted peak responses to combination of shocks to relative
consumption preferences for goods, aggregate labor supply, and the policy rate (in purple)
and peak responses in the data (in orange). Red line is model response at the posterior
mode. Black dashed lines are at steady state levels. The purple shaded areas are 68- and
95-percent credibility sets. Top left: Aggregate inflation. Orange is peak CPI inflation of
8.9 percent in June 2022. Red line is model response at posterior mode. Dashed black line is
at steady state inflation rate of 2 percent. Top right: Increase in tightness in goods sector.
V × Y ) in model
Bottom left: Peak decline in aggregate measured matching eﬀiciency (H/
and data (excluding early 2020 shutdown). Bottom right: Increase in tightness in services.

√

In other words, the shift in consumption preferences made inflation more responsive
to changes in aggregate demand, therefore increasing the risk of maintaining stimulative
In the modal case, while the model
monetary policy while tightness was already high.

38

ModelDataPercentagePoints0510Mode95%c.s.68%c.s.ModelDataV=U00.511.5Mode95%c.s.68%c.s.ModelData%dev.froms.s.-10-50Mode95%c.s.68%c.s.ModelDataV=U00.511.5Mode95%c.s.68%c.s.predicts a decline in unemployment from about 5 to about 4 percent, the monetary stimulus
needed to match inflation contributes only about a 0.15 to this reduction in unemployment.

5 Conclusion

The recent post-pandemic inflation surge highlights the potential importance of both sector-
specific shocks and labor market frictions for aggregate inflation. Sector-specific shocks
that propagate through the production network, can generate large aggregate inflationary
pressures, particularly if they hit sectors with inelastic input supply curves; such sectors
struggle to expand output and instead raise prices. I show how labor market frictions at the
sector level endogenously change the elasticity of the input supply curves faced by sectors.
Sectors with tight labor markets find it harder to hire and thus to adjust output when
additional shocks hit. As a result, labor market frictions lead to a nonlinear Phillips curve,
even if only a few sectors are constrained, making monetary policy less effective at stabilizing
output and more effective at combating inflation.

Consequently, uneven demand and supply across the economy can generate aggregate
inflation, a decline in aggregate matching eﬀiciency, and a steeper aggregate Phillips curve
qualitatively consistent with observed post-pandemic dynamics.
In addition, the model
predicts that monetary policy, while a useful stabilization tool in normal times, particularly
in the face of aggregate disturbances that affect all sectors relatively equally, may be less
effective when the economy experiences large shocks or an uneven distribution of demand
across sectors. In these circumstances, introducing sector-specific wage subsidies to ease the
transition of workers across sectors may be effective at alleviating and equalizing tightness
across the economy. The optimal policy response to sector-specific shocks in light of labor
market frictions is an interesting avenue for future work.

39

References

Abraham, K. G., & Katz, L. F. (1986). Cyclical Unemployment: Sectoral Shifts or Aggregate
Disturbances? Journal of Political Economy, 94(3), 507–522. Retrieved October 23,
2025, from https://www.jstor.org/stable/1833046

Acemoglu, D., Carvalho, V. M., Ozdaglar, A., & Tahbaz-Salehi, A. (2012). The Network
Origins of Aggregate Fluctuations. Econometrica, 80(5), 1977–2016. https://doi.org/
10.3982/ECTA9623

Afrouzi, H., Blanco, A., Drenik, A., & Hurst, E. (2024, December). A Theory of How Workers

Keep Up With Inflation. 33233. https://doi.org/10.3386/w33233

Amiti, M., Fed, N., Heise, S., Fed, N., & Karahan, F. (2023). Inflation Strikes Back: The

Role of Import Competition and the Labor Market.

Artuç, E., Chaudhuri, S., & McLaren, J. (2010). Trade Shocks and Labor Adjustment:
A Structural Empirical Approach. American Economic Review, 100(3), 1008–1045.
https://doi.org/10.1257/aer.100.3.1008

Atalay, E. (2017). How Important Are Sectoral Shocks? American Economic Journal: Macroe-

conomics, 9(4), 254–280. https://doi.org/10.1257/mac.20160353

Autor, D., Dube, A., & McGrew, A. (2023, March). The Unexpected Compression: Competi-
tion at Work in the Low Wage Labor Market. 31010. https://doi.org/10.3386/w31010
Bagga, S., Mann, L. F., Şahin, A., & Violante, G. L. (2025, May). Job Amenity Shocks and

Labor Reallocation. 33787. https://doi.org/10.3386/w33787

Baqaee, D. R. (2018). Cascading Failures in Production Networks. Econometrica, 86(5),
1819–1838. Retrieved October 25, 2022, from https://www.jstor.org/stable/44955259
Baqaee, D. R., & Farhi, E. (2019). The Macroeconomic Impact of Microeconomic Shocks:
Beyond Hulten’s Theorem. Econometrica, 87 (4), 1155–1203. https : / / doi . org / 10 .
3982/ECTA15202

Baqaee, D. R., & Farhi, E. (2022). Supply and Demand in Disaggregated Keynesian Economies
with an Application to the COVID-19 Crisis. American Economic Review, 112(5),
1397–1436. https://doi.org/10.1257/aer.20201229

Baqaee, D. R., & Rubbo, E. (2022). Micro Propagation and Macro Aggregation. 30538.

https://doi.org/10.3386/w30538

Barlevy, G., Faberman, R. J., Hobijn, B., & Şahin, A. (2023, October). Shiftin Reasons for

Beveridge-Curve Shifts.

Barnichon, R., & Shapiro, A. H. (2024). Phillips meets Beveridge. Journal of Monetary

Economics, 148, 103660. https://doi.org/10.1016/j.jmoneco.2024.103660

40

Benigno, P., & Eggertsson, G. (2023, April). It’s Baaack: The Surge in Inflation in the
2020s and the Return of the Non-Linear Phillips Curve (w31197). National Bureau
of Economic Research. Cambridge, MA. https://doi.org/10.3386/w31197

Blanchard, O., & Galí, J. (2010). Labor Markets and Monetary Policy: A New Keynesian
Model with Unemployment. American Economic Journal: Macroeconomics, 2(2), 1–
30. https://doi.org/10.1257/mac.2.2.1

Boehm, C. E., & Pandalai-Nayar, N. (2022). Convex Supply Curves. American Economic

Review, 112(12), 3941–3969. https://doi.org/10.1257/aer.20210811

Caliendo, L., Dvorkin, M., & Parro, F. (2019). Trade and Labor Market Dynamics: General
Equilibrium Analysis of the China Trade Shock. Econometrica, 87 (3), 741–835. https:
//doi.org/10.3982/ECTA13758

Cardoza, M., Grigoli, F., Pierri, N., & Ruane, C. (2022). Worker Mobility in Production

Networks.

Carvalho, C., Nechio, F., & Tristão, T. (2021). Taylor rule estimation by OLS. Journal of
Monetary Economics, 124, 140–154. https://doi.org/10.1016/j.jmoneco.2021.10.010
Carvalho, V. M., & Tahbaz-Salehi, A. (2019). Production Networks: A Primer. Annual Re-
view of Economics, 11(1), 635–663. https://doi.org/10.1146/annurev- economics-
080218-030212

Chetty, R. (2012). Bounds on Elasticities With Optimization Frictions: A Synthesis of Micro
and Macro Evidence on Labor Supply. Econometrica, 80(3), 969–1018. https://doi.
org/10.3982/ECTA9043

Comin, D. A., Johnson, R. C., & Jones, C. J. (2023, April). Supply Chain Constraints and

Inflation. 31179. https://doi.org/10.3386/w31179

Di Giovanni, J., Kalemli-Özcan, Ş., Silva, A., & Yıldırım, M. A. (2023). Quantifying the
Inflationary Impact of Fiscal Stimulus under Supply Constraints. AEA Papers and
Proceedings, 113, 76–80. https://doi.org/10.1257/pandp.20231028

Diamond, P. A. (1982a). Aggregate Demand Management in Search Equilibrium. Journal of

Political Economy.

Diamond, P. A. (1982b). Wage Determination and Eﬀiciency in Search Equilibrium. The
Review of Economic Studies, 49(2), 217–227. https://doi.org/10.2307/2297271
di Giovanni, J., Kalemli-Özcan, Ṣ., Silva, A., & Yildirim, M. A. (2023, November). Pandemic-
Era Inflation Drivers and Global Spillovers. 31887. https://doi.org/10.3386/w31887
Faccini, R., & Melosi, L. (2025). Job-to-Job Mobility and Inflation. The Review of Economics

and Statistics, 107 (4), 1027–1041. https://doi.org/10.1162/rest_a_01312

41

Ferrante, F., Graves, S., & Iacoviello, M. (2023). The inflationary effects of sectoral reallo-
cation. Journal of Monetary Economics, 140, S64–S81. https://doi.org/10.1016/j.
jmoneco.2023.03.003

Fujita, S., Moscarini, G., & Postel-Vinay, F. (2024). Measuring Employer-to-Employer Re-
allocation. American Economic Journal: Macroeconomics, 16(3), 1–51. https://doi.
org/10.1257/mac.20210076

Gertler, M., & Trigari, A. (2009). Unemployment Fluctuations with Staggered Nash Wage
Bargaining. Journal of Political Economy, 117 (1), 38–86. https://doi.org/10.1086/
597302

Gitti, G. (2024, January). Nonlinearities in the Regional Phillips Curve with Labor Market

Tightness.

Guerrieri, V., Lorenzoni, G., Straub, L., & Werning, I. (2022). Macroeconomic Implications
of COVID-19: Can Negative Supply Shocks Cause Demand Shortages? American
Economic Review, 112(5), 1437–1474. https://doi.org/10.1257/aer.20201063
Hall, R. (2018, May). New Evidence on the Markup of Prices over Marginal Costs and the
Role of Mega-Firms in the US Economy (w24574). National Bureau of Economic
Research. Cambridge, MA. https://doi.org/10.3386/w24574

Hall, R. E. (2005). Employment Fluctuations with Equilibrium Wage Stickiness. The Amer-
ican Economic Review, 95(1), 50–65. Retrieved May 11, 2022, from http://www.
proquest.com/docview/233027942/abstract/10968E318BD94539PQ/1

Hazell, J., Herreño, J., Nakamura, E., & Steinsson, J. (2022). The Slope of the Phillips Curve:
Evidence from U.S. States*. The Quarterly Journal of Economics, 137 (3), 1299–1344.
https://doi.org/10.1093/qje/qjac010

Hulten, C. R. (1978). Growth Accounting with Intermediate Inputs. The Review of Economic

Studies, 45(3), 511–518. https://doi.org/10.2307/2297252
Humlum, A. (2021). Robot Adoption and Labor Market Dynamics, 69.
Jones, C. I. (2011). Misallocation, Economic Growth, and Input-Output Economics. 16742.

https://doi.org/10.3386/w16742

La’O, J., & Tahbaz-Salehi, A. (2022). Optimal Monetary Policy in Production Networks.

Econometrica, 90(3), 1295–1336. https://doi.org/10.3982/ECTA18627

Lilien, D. M. (1982). Sectoral Shifts and Cyclical Unemployment. Journal of Political Econ-
omy, 90(4), 777–793. Retrieved October 23, 2025, from https : / / www . jstor . org /
stable/1831352

Lorenzoni, G., & Werning, I. (2024). Discussion of Revisiting the Phillips and Beveridge

Curves: Insights from the 2020s Inflation Surge. Jackson Hole Symposium.

42

McLeay, M., & Tenreyro, S. (2020). Optimal Inflation and the Identification of the Phillips

Curve. NBER Macroeconomics Annual.

Michaillat, P., & Saez, E. (2024, October 27). Beveridgean Phillips Curve. arXiv: 2401.12475

[econ]. https://doi.org/10.48550/arXiv.2401.12475

Minton, R., & Wheaton, B. (2023). Delayed Inflation in Supply Chains: Theory and Evidence.

SSRN Electronic Journal. https://doi.org/10.2139/ssrn.4470302

Mortensen, D. T. (1982a). The Matching Process as a Noncooperative Bargaining Game. In
The Economics of Information and Uncertainty (pp. 233–258). University of Chicago
Press. Retrieved May 6, 2022, from https://www.nber.org/books- and- chapters/
economics-information-and-uncertainty/matching-process-noncooperative-bargaining-
game

Mortensen, D. T. (1982b). Property Rights and Eﬀiciency in Mating, Racing, and Related
Games. American Economic Review, 72(5), 968. Retrieved November 1, 2023, from
https://search.ebscohost.com/login.aspx?direct=true&db=bth&AN=4496919&
site=ehost-live&scope=site&authtype=ip,sso&custid=rock

Moscarini, G., & Postel-Vinay, F. (2023). The Job Ladder: Inflation vs. Reallocation.
Pasten, E., Schoenle, R., & Weber, M. (2020). The propagation of monetary policy shocks
in a heterogeneous production economy. Journal of Monetary Economics, 116, 1–22.
https://doi.org/10.1016/j.jmoneco.2019.10.001

Pissarides, C. A. (1984). Search Intensity, Job Advertising, and Eﬀiciency. Journal of Labor

Economics, 2(1), 128–143. https://doi.org/10.1086/298026

Pissarides, C. A. (1985). Short-Run Equilibrium Dynamics of Unemployment, Vacancies,
and Real Wages. American Economic Review, 75(4), 676–690. Retrieved November
1, 2023, from https://search.ebscohost.com/login.aspx?direct=true&db=bth&AN=
4500845&site=ehost-live&scope=site&authtype=ip,sso&custid=rock

Rubbo, E. (2023). Networks, Phillips Curves, and Monetary Policy. Econometrica, 91(4),

1417–1455. https://doi.org/10.3982/ECTA18654

Rubbo, E. (2024). What drives inflation? Lessons from disaggregated price data.
Şahin, A., Song, J., Topa, G., & Violante, G. L. (2014). Mismatch Unemployment. American
Economic Review, 104(11), 3529–3564. https://doi.org/10.1257/aer.104.11.3529

Schüle, F., & Sheng, H. (2024). Unemployment in a Production Network (Working Paper).
Shimer, R. (2005). The Cyclical Behavior of Equilibrium Unemployment and Vacancies.
The American Economic Review, 95(1), 25–49. Retrieved May 10, 2022, from https:
//www.proquest.com/docview/233026649/abstract/92577F8A4D754A9FPQ/1

43

44

A All Equilibrium Conditions

The following conditions hold for each sector i ∈ J:

JX

ϵy −1
ϵy
i,t +

1
ϵy

ω

ij X

!! ϵy
ϵy −1

ϵy −1
ϵy
ij,t

1
ϵy

ω

in N

Yi,t = Ai,t

= µi,t (qi,t − ri,t) + Et

(cid:2)

i=1

SDFt|t+1(1 − si,t+1)µi,t+1ri,t

(cid:3)

ri,t

Wi,t
Pt
qi,t
ri,t

µi,t

Pj,tXij,t
M Ci,tYi,t

=

M Ci,t
Pt

1
ϵy

in A
β
(cid:18)

ϵy −1
ϵy
i,t N

− 1
ϵy
i,t Y

1
ϵy
i,t
(cid:19) ϵy −1
ϵy

= (βixωij)

(cid:18)

Xij,t
Yi,t

1
ϵy

Ai,t
(cid:19)−ϵd

C agg
t

∀ j ∈ J

(cid:18)

(cid:19)

Πi,t
Π

− 1

Ci,t = αi,t

Πi,t
Π

=

ϵ
ψp

(cid:18)

Pi,t
Pt
M Ci,t
Pt

− ϵ − 1
ϵ

Pi,t
Pt

(cid:19)

(cid:20)

(cid:18)

+ Et

SDFt|t+1

(cid:19)

Πi,t+1
Π

− 1

Πi,t+1
Π

Yi,t+1
Yi,t

Li,t/Lt
Li,t−1/Lt−1

(cid:21)

− 1

+

(cid:20)

1
2

!

(cid:21)2

Li,t/Lt
Li,t−1/Lt−1

#

− 1

"

(cid:20)

Zt+1

Li,t+1/Lt+1
Li,t/Lt

(cid:21)

− 1
(cid:19)

L2

i,t+1/Lt+1
L2
i,t/Lt

+ βψL,iEt
(cid:18)

− Ξi,tLφ
i,t

+ βEt [(1 − fi,t+1)(1 − si,t+1)Ξi,t+1]

(47)

#

SDFt|t+1 = βEt

Πi,t = Πagg

t

"(cid:18)

(cid:18)

(cid:19)−σ Zt+1
Zt
(cid:19)

C agg
t+1
C agg
t
Pi,t/Pt
Pi,t−1/Pt−1
JX

Yi,t = Ci,t +

Xji,t

j=1

 (cid:20)

δt = Ξi,tfi,t − ψL,iZt

−σ Wi,t
Pt

Ξi,t = Zt

(C agg
t

)

JX

Lt =

Li,t

i=1

Li,t = Ni,t + ri,tVi,t
Li,t = (1 − si,t)Li,t−1 + Ui,t

qi,t =

fi,t =

Hi,t
Vi,t
Hi,t
Ui,t

45

(36)

(37)

(38)

(39)

(40)
(cid:21)

(41)

(42)

(43)

(44)

(45)

(46)

(48)

(49)

(50)

(51)

(52)

(cid:1)

(53)

(54)

(55)

(56)

(57)

(58)

(59)

(60)

(61)

(62)

(63)

(64)

θi,t =

Vi,t
Ui,t
(cid:0)

(cid:1)− 1
ηi

U

−ηi
i,t

Hi,t = ζi,t

−ηi
i,t + V
Ni,t + ri,tVi,t = Hi,t + (1 − si,t)(Ni,t−1 + ri,t−1Vi,t−1)
1
−σ − κ ri,t
)
qi,t

Ztχi,tLφ
i,t

Zt (C agg

Wi,t
Pt

=

(cid:0)

t

− βEt [(1 − fi,t+1)(1 − si,t+1)Ξi,t+1]

! ϵd
−1
ϵd

ϵd

−1
ϵd

i,t

Aggregate variables satisfy:

C agg

t =

JX

1
ϵd
i,t C

α

H agg

t =

V agg
t =

U agg

t =

θagg
t =

i=1

JX

Hi,t

i=1
JX

i=1
JX

Vi,t

Ui,t

i=1
V agg
t
U agg
t

1 = β(1 + it)Et

"(cid:18)

(cid:16)

(cid:19)−σ Zt+1
Zt

C agg
t+1
C agg
t

1
Πagg
t
(cid:17)
1−ρi

#

Mt

(1 + it) = (1 + it−1)ρi

[Πagg
t

]ϕπ [Y agg
t

]ϕy

t = C agg
Y agg

t

46

Finally, the shock processes are given by:

si,t = sρs
Ai,t = AρA
αi,t = αρα

i,t−1εs,i,t
i,t−1εA,i,t
i,t−1εα,i,t ∀ i ∈ J − 1

1
ϵd

J,t = 1 −

α

J−1X

1
ϵd
i,t

α

ri,t = rρr
ζi,t = ζ ρζ
χi,t = χρχ
Zt = Z ρZ
Mt = M ρM

i=1
i,t−1εr,i,t
i,t−1εζ,i,t
i,t−1εχ,i,t
t−1εZ,t
t−1εm,t

(65)

(66)

(67)

(68)

(69)

(70)

(71)

(72)

(73)

A.1 Adding vacancy Adjustment Costs and Price Indexation to

Past Inflation

L (·) = Et

(cid:26)

∞X

(cid:20)

SDFt|t+s

(cid:19)−ϵ

(cid:18)

Pi,t+s(z)
Pi,t+s

Pi,t+s(z)
Pt+s
(cid:18)

Xij,t+s(z) − ψp
2

Π1−ρπ Πρπ

Pi,t+s(z)
i,t−1Pi,t+s−1(z)
#

X

−

j

s=0
Pj,t
Pt
(cid:20)

+ λi,t+s

Ai,t+s
(cid:20)

1
ϵy
ix

β

"

X

j

1
ϵy

ω

ij Xij,t+s(z)

ϵy −1
ϵy

1
ϵy

+ β

in Ni,t+s(z)

Yi,t+s − Wi,t+s
Pt+s
(cid:19)
2

− 1

Yi,t+s

[Ni,t+s(z) + ri,t+sVi,t+s(z)]

(cid:21)

! ϵy
ϵy −1

(cid:18)

−

ϵy −1
ϵy

Pi,t+s(z)
Pi,t+s

(cid:19)−ϵ

(cid:21)

Yi,t+s

(cid:21)(cid:27)

+ µi,t+s

(1 − si,t+s) (Ni,t+s−1(z) + ri,t+s−1Vi,t+s−1(z)) + (qi,t+s − ri,t+s) Vi,t+s(z) − Ni,t+s(z)

FOCs

∂L
∂Pi,t(z)

(cid:19)−ϵ

(cid:18)

1
Pt

Pi,t(z)
Pi,t

= (1 − ϵ)
( (cid:18)

Yi,t + ϵλi,t
(cid:19)

− ψp

(cid:20)

Pi,t(z)
Π1−ρπ Πρπ
i,t−1Pi,t−1(z)
(cid:18)

− 1

(cid:18)

(cid:19)−ϵ Yi,t
Pi,t(z)

Pi,t(z)
Pi,t

1

Π1−ρπ Πρπ

i,t−1Pi,t−1(z)

Yi,t

− Et

SDFt|t+s

Pi,t+1(z)

Π1−ρπ Πρπ

i,t Pi,t(z)

Pi,t+1(z)

Π1−ρπ Πρπ

i,t Pi,t(z)2 Yi,t+1

(cid:19)

− 1

(cid:21) )

= 0

Now I assume a symmetric equilibrium, where we can drop the z indexation. Starting

47

with (74), this implies
(cid:18)

(cid:19)

Πi,t
Π1−ρπ Πρπ

i,t−1

− 1

1
Π1−ρπ Πρπ
i,t−1Pi,t−1

Yi,t =

(cid:20)

1
ψp
(cid:20)

(1 − ϵ)

1
Yi,t + ϵλi,t
Pt
(cid:18)

(cid:21)

Yi,t
Pi,t

+ Et

SDFt|t+s

Πi,t+1
Π1−ρπ Πρπ

i,t Pi,t

− 1

(cid:19)

(cid:21)

Pi,t+1
Π1−ρπ Πρπ

i,t P 2
i,t

Yi,t+1

Dividing by Yi,t and multiplying by Pi,t gives
(cid:20)
(cid:19)

(cid:18)

Πi,t
Π1−ρπ Πρπ

i,t−1

− 1

Πi,t
Π1−ρπ Πρπ

i,t−1

=

ϵ
ψp
(cid:20)

λi,t − ϵ − 1
ϵ
(cid:18)

+ Et

SDFt|t+s

(cid:21)

Pi,t
Pt
Πi,t+1
Π1−ρπ Πρπ
i,t

(cid:19)

− 1

(cid:21)

Πi,t+1
Π1−ρπ Πρπ
i,t

Yi,t+1
Yi,t+1

I add vacancy adjustment costs by modeling ri,t as a quadratic function that depends on

the growth of firm level vacancy postings from one period to the next.

B First Order Approximation to the Phillips Curve

max EtSDFt+s|t

Di,t+s(z)
Pt+s

Di,t(z) =

Pi,t(z)
Pt

(cid:18)

Pi,t(z)
Pi,t

(cid:19)−ϵ

Yi,t − Wi,t
Pt

− Wi,t
Pt

(cid:2)

Ni,t(z) + N r

i,t(z)

(cid:3)

−

JX

j=1

Pj,t
Pt

Xij,t(z) − ψp
2

(cid:18)

Pi,t(z)
ΠPi,t−1(z)

(cid:19)

2

− 1

Yi,t

! ϵy
ϵy −1

Ni,t(z) + ri,tVi,t(z) = qi,tVi,t(z) + (1 − si,t)(Ni,t−1 + ri,t−1Vi,t−1)

Yi,t(z) = Ai,t

1
ϵy

ω

in Ni,t(z)

1
ϵy

ω

ij Xij,t(z)

ϵy −1
ϵy

ϵy −1
ϵy +

JX

j=1

48

Lagrangian:

L (·) = Et

(cid:26)

∞X

(cid:20)

SDFt|t+s

Pi,t+s(z)
Pt+s
(cid:18)

(cid:18)

Pi,t+s(z)
Pi,t+s

(cid:19)−ϵ

Yi,t+s − Wi,t+s
Pt+s
(cid:21)
(cid:19)

2

[Ni,t+s(z) + ri,t+sVi,t+s(z)]

X

−

j

s=0
Pj,t
Pt
(cid:20)

+ λi,t+s

Ai,t+s
(cid:20)

j

Xij,t+s(z) − ψp
2
 "

X

Pi,t+s(z)
ΠPi,t+s−1(z)
#

− 1

Yi,t+s

1
ϵy

ω

ij Xij,t+s(z)

ϵy −1
ϵy

1
ϵy

+ ω

in Ni,t+s(z)

! ϵy
ϵy −1

(cid:18)

−

ϵy −1
ϵy

Pi,t+s(z)
Pi,t+s

(cid:19)−ϵ

(cid:21)

Yi,t+s
(cid:21)(cid:27)

+ µi,t+s

(1 − si,t+s) (Ni,t+s−1 + ri,t+s−1Vi,t+s−1) + (qi,t+s − ri,t+s) Vi,t+s(z) − Ni,t+s(z)

FOCs

∂L
∂Pi,t(z)

= (1 − ϵ)
(cid:26)(cid:18)

− ψp

(cid:18)

1
Pt

Pi,t(z)
Pi,t

(cid:19)−ϵ

(cid:19)

Pi,t(z)
ΠPi,t−1(z)

− 1

(cid:18)

Yi,t + ϵλi,t

Pi,t(z)
Pi,t

(cid:19)−ϵ Yi,t
Pi,t(z)

(cid:20)

(cid:18)

1
ΠPi,t−1(z)
X

(  "

Yi,t − Et

SDFt|t+s
#

(cid:19)

Pi,t+1(z)
ΠPi,t(z)

− 1

(74)

Pi,t+1(z)
ΠPi,t(z)2 Yi,t+1
! ϵy
ϵy −1

−1

(cid:21)(cid:27)

= 0

∂L
∂Xij,t(z)

= − Pj,t
Pt

+ λi,tAi,t

ϵy
ϵy − 1

1
ϵy

ω

ij Xij,t+s(z)

ϵy −1
ϵy

1
ϵy

+ ω

in Ni,t+s(z)

ϵy −1
ϵy

(75)

! ϵy
ϵy −1

−1

ϵy −1
ϵy

(76)

(77)

(78)

j

)

1
× ω
ϵy
ij

ϵy − 1
ϵy

Xij,t(z)

ϵy −1
ϵy

−1

= 0

∂L
∂Ni,t(z)

= − Wi,t
Pt

+ λi,tAi,t

ϵy
ϵy − 1

(  "

X

j

)

1
ϵy

ω

ij Xij,t+s(z)

#

ϵy −1
ϵy

1
ϵy

+ ω

in Ni,t+s(z)

Nij,t(z)

ϵy −1
ϵy

−1

− µi,t = 0

ri,t + (qi,t − ri,t) µi,t = 0

1
ϵy
in

× ω

ϵy − 1
ϵy
= − Wi,t
Pt

∂L
∂Vi,t(z)

Implies,

µi,t =

ri,t
qi,t − ri,t

µi,t = λi,tω

1
ϵy

in A

Wi,t
Pt
ϵy −1
ϵy
i,t Y

1
ϵy

i,t N

− 1
ϵy

i,t

− Wi,t
Pt

49

Combining implies,

⇒ µi,t +

qi,t − ri,t
ri,t
⇒ qi,t
ri,t

µi,t = λi,tω

1
ϵy

in A

µi,t = λi,tω

1
ϵy

in A

µi,t = λi,tω

1
ϵy

in A

− qi,t − ri,t
ri,t

µi,t

ϵy −1
ϵy
i,t Y

1
ϵy

i,t N

− 1
ϵy

i,t

ϵy −1
ϵy
i,t Y

1
ϵy

i,t N

− 1
ϵy

i,t

ϵy −1
ϵy
i,t Y

1
ϵy

i,t N

− 1
ϵy

i,t

Or alternatively,

Similarly,

qi,t
qi,t − ri,t

Wi,t
Pt

1
ϵy

= λi,tω

in A

ϵy −1
ϵy
i,t Y

1
ϵy

i,t N

− 1
ϵy

i,t

Pj,t
Pt

= (ωij)

1

ϵy λi,tA

ϵy −1
ϵy
i,t Y

1
ϵy

i,t X

− 1
ϵy
ij,t

Combining, we can write all Xij,t and Ni,t in terms of Xii,t:

qi,t
qi,t − ri,t

Wi,t
Pt

=

Pi,t
Pt

⇒ N

1
ϵy

i,t =

⇒ Ni,t =

Pi,t
Wi,t
(cid:18)
ωin
ωii

(cid:18)

(cid:19) 1
ϵy

ωin
ωii
qi,t − ri,t
qi,t
(cid:19) (cid:18)

Pi,t
Wi,t

(cid:19) 1
ϵy

(cid:18)

Xii,t
Ni,t

(cid:18)

(cid:19) 1
ϵy

ωin
ωii
qi,t − ri,t
qi,t

1
ϵy
ii,t

X
(cid:19)

ϵy

Xii,t

And similarly,

Xij,t =

(cid:18)

ωij
ωii

(cid:19) (cid:18)

(cid:19)

ϵy

Pi,t
Pj,t

Xii,t

Now plugging into the production function for optimal input choices,

(cid:18)(cid:18)

0

@ω

1
ϵy
in

ωin
ωii

(cid:19) (cid:18)

Yi,t = Ai,t

Pi,t
Wi,t

qi,t − ri,t
qi,t

(cid:19)

ϵy

Xii,t

(cid:19) ϵy −1
ϵy

+

JX

j=1

1
ϵy
ij

ω

(cid:18)(cid:18)

(cid:19) (cid:18)

(cid:19)

ϵy

Pi,t
Pj,t

ωij
ωii

Xii,t

(cid:19) ϵy −1
ϵy

ϵy
ϵy −1

1

A

(cid:18)

= Ai,t

ωin

Wi,t

qi,t
qi,t − ri,t

(cid:19)

1−ϵy

+

JX

j=1

ωijP 1−ϵy

j,t

! ϵy

ϵy −1 1
ωii

P ϵy

i,t Xii,t

50

We can rearrange this to

Xii,t = ωiiA−1

i,t

(cid:19)−ϵy

(cid:18)

Pi,t
Θi,t

Yi,t

(cid:18)

(cid:16)

Where, Θi,t =

ωin

Wi,t

(cid:17)

1−ϵy

P

+

qi,t
qi,t−ri,t

J

j=1 ωijP 1−ϵy

j,t

(cid:19) 1
1−ϵy

. This lets us write

Ni,t = ωinA−1

i,t

(cid:18)

Xij,t = ωijA−1

i,t

We can then write the cost in terms of Yi,t,

(cid:18)

Costi,t = Wi,t

qi,t
qi,t − ri,t
(cid:18)

Ni,t − ri,t
qi,t − ri,t
(cid:19)

1−ϵy

qi,t
qi,t − ri,t

= ωinA−1

i,t

Wi,t

The marginal cost is then,

(cid:19)−ϵy

Yi,t

qi,t
qi,t − ri,t
(cid:19)−ϵy

(cid:18)

Wi,t
Θi,t
Pj,t
Θi,t

Yi,t

(cid:19)

JX

(Ni,t−1 + ri,t−1Vi,t−1)

+

Pj,tXij,t

j=1

Θϵy

i,tYi,t +

JX

j=1

ωijA−1

i,t P 1−ϵy

j,t Θϵy

i,tYi,t + ...

M Ci,t = A−1

i,t Θϵy

i,t

(cid:18)

ωin

Wi,t

qi,t
qi,t − ri,t

(cid:19)

1−ϵy

+

JX

j=1

ωijP 1−ϵy

j,t

!

= A−1

i,t Θϵy

i,tΘ1−ϵy
i,t = A−1
(cid:18)

i,t Θi,t

= A−1
i,t

ωin

Wi,t

qi,t
qi,t − ri,t

(cid:19)

1−ϵy

+

JX

j=1

ωijP 1−ϵy

j,t

! 1
1−ϵy

All other conditions are the same as in the full model.

On the household side, assume there are no labor rellocation costs, so that the household

51

problem is

"

∞X

max Et

Zt+s

s=0

Ct+s(z)1−σ
1 − σ

−

(cid:18)

JX

i=1

χi,t+s

Li,t+s(z)1+φ
1 + φ

(cid:19) #

JX

s.t. Pt+sCt+s(z) + Bt+s(z) = (1 + it+s)Bt+s−1(z) +

Wi,t+sLi,t+s(z) + Tt+s(z)

Li,t+s(z) = (1 − fi,t+s)(1 − si,t+s)Li,t+s−1 + fi,t+sLi,t+s(z)

i=1

and,

JX

i=1

Li,t+s(z) ≤ Lt+s

The Lagrangian is

L = Et

∞X

s=0
(cid:20)

(cid:18)

"

βs

Zt+s

Ct+s(z)1−σ
1 − σ

−

(cid:18)

JX

i=1

χi,t+s

Li,t+s(z)1+φ
1 + φ

(cid:19)#

+ υt+s

(1 + it+s−1)Bt+s−1(z) +

JX

i=1

Wi,t+sLi,t+s(z) + Tt+s(z) − Pt+sCt+s(z) − Bt+s(z)

(cid:21)

Ξi,t+s [(1 − fi,t+s)(1 − si,t+s)Li,t+s−1 + fi,t+sLi,t+s(z) − Li,t+s(z)]
"

#

+

JX

i=1

+ δt+s

Lt+s −

JX

Li,t+s

i=1

And the first order conditions are

∂L
∂Ct+s(z)
∂L
∂Bt+s(z)
∂L
∂Li,t+s(z)
∂L
∂Li,t+s(z)

Which implies,

= βsEt
(cid:2)

= Et

(cid:20)

= Et

= Et

= 0

(cid:3)

= 0

(cid:2)

Zt+sCt+s(z)−σ − υt+sPt+s

(cid:3)

βs+1υt+s+1(1 + it+s) − βsυt+s

(cid:21)

= 0

βsΞi,t+sfi,t+s − βsδt+s
(cid:2)

βsυt+sWi,t+s − βsZtχi,t+sLφ

i,t+s

− βsΞi,t+s

(cid:3)

= 0

(cid:18)

Zt

C −σ
t

Wi,t
Pt

Ξi,tfi,t = δt
(cid:19)

− χi,tLφ
i,t

= Ξi,t

52

The Nash bargaining solution is then

(cid:18)

(cid:19)

Wi,t
Pt

− χi,tLφ
i,t

Wi,t =

Pt
ZtC −σ
Wi,t = Wi,t − χi,tLφ

Zt

t

C −σ
t

i,tC σ

t Pt

(cid:18)

⇒

1 − κ

κ

κ

ri,t
qi,t
ri,t
qi,t
(cid:19)

ri,t
qi,t
⇒ Wi,t
Pt

Wi,t = χi,tLφ

i,tC σ

t Pt

=

1
1 − κ ri,t
qi,t

χi,tLφ

i,tC σ

t

Suppose

"

"

1
1 − κ ri,t
qi,t

1
1 − κ ri,t
qi,t

Wi,t =

⇒ Wi,t
Pt

=

#

1−ρw

χi,tLφ

i,tC σ

t Pt
#

1−ρw (cid:18)

χi,tLφ

i,tC σ

t

W ρw

i,t−1

(cid:19)

ρw

Wi,t−1
Pt−1

Π−1
t

B.1 First order approx

To first order, the NKPC is

πi,t =

ϵ − 1
ψp,i

(mci,t − pi,t) + βEtπi,t+1

For a general matching function m(Ui, Vi), the job-finding rate is

qi,t = mi,t (1, θi,t)

(θi,t − θi)

qi + (qi,t − qi) = qi +
qi,t − qi
qi

=

∂qi
∂θi
qi,t = E qi
θi

∂qi
∂θi
θi,t − θi
θi
θi
qi
bθi,t

Where E qi
θi

is the elasticity of the vacancy-filling rate to changes in tightness.

Now for the marginal cost,

53

!

(cid:19)

1−ϵy

Pj,t
Pt
(cid:19)

1−ϵy

Wi,t
Pt
(cid:19)−ϵy

qi,t
qi,t − ri,t
(cid:18)

(cid:19)

1−ϵy

JX

(cid:18)

+

ωij

(cid:18)

M Ci,t
Pt
Wi,t
Pt

− M Ci
P
(cid:19)
− Wi
P

(cid:19)

j=1

(cid:18)

=

M Ci
Pt

− (1 − ϵy)

(cid:18)

(cid:18)

M Ci,t
Pt

M Ci
Pt

(cid:19)

1−ϵy

= A

−(1−ϵy)
i,t

ωin

(cid:18)

(cid:19)

1−ϵy

(cid:18)

+ (1 − ϵy)
(cid:18)

M Ci
Pt
(cid:19)

+ (1 − ϵy)

1
Ai

ωin
(cid:18)

Wi
P

qi
qi − ri
(cid:19)

+ (1 − ϵy)

1
Ai
(cid:18)

Wi
P

qi
qi − ri
(cid:19)
1−ϵy

1−ϵy P
Wi

1−ϵy

qi

⇒ (1 − ϵy)

+ (1 − ϵy)

+ (1 − ϵy)

1
Ai

1
Ai

M Ci
Pt

(cid:18)

ωin

ωin

(cid:18)

(cid:19)

1−ϵy

(cid:19)

1−ϵy

Wi
P

qi
qi − ri

Wi
P

qi
qi − ri

mci,t = −(1 − ϵy)

(qi − ri)2 (ri,t − ri) + (1 − ϵy)
M Ci
Pt

1−ϵy

(cid:19)

(cid:18)

ai,t
(cid:18)

wi,t − (1 − ϵy)

1
Ai

ωin

Wi
P

riqi

(qi − ri)2 ri,t + (1 − ϵy)

1
Ai

(cid:19)

1−ϵy

(cid:18)

Pj
P

qi
qi − ri
JX

ωij

j=1

riqi
(qi − ri)2 qi,t
(cid:19)

1−ϵy

pj,t

(cid:18)

(cid:19)

M Ci
Pt
1−ϵy

(cid:19)

− (1 − ϵy)
(cid:18)

Wi
P
(cid:18)

qi
qi − ri
(cid:19)

(cid:18)

Pj
P

1−ϵy P
Pj

ri

(Ai,t − Ai)

1−ϵy 1
Ai
(qi − ri)2 (qi,t − qi)
(cid:19)
Pj,t
Pt

− Pj
P

ωin

1
Ai
JX

1
Ai

ωij

j=1

From the first order conditions,

(cid:18)

(cid:18)

Ni
Yi

(cid:19) ϵy −1
ϵy

= Ωin

(cid:19)− 1

ϵy

WiNi
M CiYi

ϵy −1
ϵy

= ω

in A
i

1
ϵy

1
ϵy

qi
qi − ri
⇒ qi

qi − ri

ϵy −1
ϵy

Ni
Yi
(1−ϵy )(ϵy −1)
ϵy

= ω

in A
i

1−ϵy
ϵy
in A
i

= ω

Wi
M Ci
(cid:21)
1−ϵy

(cid:21)

1−ϵy

(cid:18)

(cid:19) ϵy −1
ϵy

Ni
Yi
(cid:18)
(1−ϵy )(ϵy −1)
ϵy

(cid:19) ϵy −1
ϵy

=

ωin
Ai

1−ϵy
ϵy
in A
i

ω

= A

−ϵy
i ω

1
ϵy

in A
i

ϵy −1
ϵy

Ni
Yi
(cid:19) ϵy −1
ϵy

(cid:18)

Ni
Yi

(cid:20)

(cid:20)

⇒

⇒ ωin
Ai

qi
qi − ri

Wi
M Ci

qi
qi − ri

Wi
M Ci

= A

−ϵy
i Ωin = Ωin

So we can rewrite the marginal cost as

mci,t = −ai,t + Ωin

(cid:20)

wi,t − ri

qi − ri

qi
qi − ri

(cid:21)

(qi,t − ri,t)

+

JX

j=1

Ωijpj,t

54

Or in terms of tightness,

mci,t = −ai,t + Ωin

(cid:20)

wi,t − ri

qi − ri

qi
qi − ri

(cid:21)

(E qi
θi

bθi,t − ri,t)

+

JX

j=1

Ωijpj,t

Stacking over sectors, and assuming wages are fully rigid,

(η bθt + rt) + Ωxpt
mct = −at + Ωn r(Q − r)−1Q(Q − r)−1
}

|

{z
ΓQ

Where η is a diagonal matrix with the negative of the elasticity of the vacancy-filling rate
to changes in tightness on the diagonal. Stacking the sector level PC over sectors gives

πt = λ (mct − pt) + βEtπt+1

Where λ is a diagonal matrix capturing pricing frictions in each sector. Plugging into for
the marginal cost gives

πt = λ (Ωx − I) pt + ΩnΓQη bθt + βEtπt+1

+ ΩnΓQrt − at

Finally, using pt = πt + pt−1 − 1πagg

t

, we can rewrite the sector level NKPC as

πt = λ (Ωx − I) (πt − 1πagg

t

) + ΩnΓQη bθt + βEtπt+1

+ ΩnΓQrt − at + λ (Ωx − I) pt−1

Now from the the aggregate consumption price index,

Pt
Pt−1

=

Π1−ϵd
t

=

(cid:19)

1−ϵd

! 1

1−ϵd

Pi,t
Pi,t−1

Pi,t−1
Pt−1
(cid:19)

Pi,t−1
Pt−1

1−ϵd

JX

(cid:18)

αi,t

(cid:18)

i=1

JX

αi,t

Πi,t−1

i=1

55

To first order,

Π1−ϵd + (1 − ϵd)Π1−ϵd

Πt − Π
Π

= Π1−ϵd +

JX

i=1

+

JX

i=1

Which implies,

(cid:18)

(cid:19)

Π1−ϵdαi

(cid:18)

Pi
P
(cid:19)

1−ϵd αi,t − αi

αi

(1 − ϵd)Π1−ϵdαi

(cid:18)

1−ϵd

Πi,t − Π
Π

+

Pi,t/Pt − Pi/P
Pi/P

(cid:19)

Pi
P

πagg
t =

(cid:20)

JX

Ωd,i
(cid:20)

i=1

1
1 − ϵd

(cid:21)

bαi,t + πi,t + pi,t−1

(cid:21)

1
1 − ϵd

αt

= Ω′
d

πt + pt−1 +

Plugging in to the sectoral PC gives

πt = λ (Ωx − I) (I − 1Ω′

+ ΩnΓQrt − at + λ (Ωx − I) (I − 1Ω′

λ (Ωx − I) 1Ω′

dαt

d) πt + ΩnΓQη bθt + βEtπt+1
d) pt−1 − 1
−1 [ΩnΓQηθt + βEtπt+1] + υsec

1 − ϵd

d)]

t

⇒ πt = [I − λ (Ωx − I) (I − 1Ω′

Where

t = [I − λ (Ωx − I) (I − 1Ω′
υsec

d)]

"

−1

ΩnΓQrt − at + λ (Ωx − I) (I − 1Ω′

d) pt−1 − 1

1 − ϵd

#

λ (Ωx − I) 1Ω′

dαt

Finally, using the expression for aggregate inflation, the aggregate NKPC is

πagg
t = Γθθt + ΓπEtπt+1 + υt

56

Where,

Γθ = Ω′
Γπ = βΩ′
υt = Ω′
+ Ω′
+ Ω′

d [I − λ (Ωx − I) (I − 1Ω′
d)]
d [I − λ (Ωx − I) (I − 1Ω′
d [I − λ (Ωx − I) (I − 1Ω′
d)]
d [I − λ (Ωx − I) (I − 1Ω′
d)]
d [I − λ (Ωx − I) (I − 1Ω′
d)]

−1 ΩnΓQη
−1
d)]
−1 [I + λ (Ωx − I) (I − 1Ω′
−1 [ΩnΓQrt − at]
−1 [I − λ (Ωx − I) 1Ω′

d] αt

d)] pt−1

C Unemployment Based Measure of Tightness

In this section, I show that the patterns I highlight in section 4.1 hold when using the more
conventional measure of tightness V
U .

57

(a) Real PCE and PCE Inflation.

(b) labor market tightness and PCE Inflation.

Figure 12: Left Axis (both subplots): Dashed lines are changes in real PCE relative to January 2020, by
major expenditure category: goods (blue) and services (orange). Right Axis: Solid lines are year-over-year
changes in the PCE price index, by major category: goods (blue) and services (orange). The gray shaded
area indicates the period from the start of the inflation surge in early 2021, to the first Federal Reserve rate
hike in March 2022.

58

201920202021202220232024202595100105110115RealConsumption(Jan2020=100)PCEInﬂationandRealConsumption(2019-Present)GoodsRealPCEServicesRealPCE−20246810PCEInﬂation(YoY%)GoodsInﬂationServicesInﬂation201920202021202220232024202550100150200250300Tightness(Jan2020=100)NormalizedTightnessandInﬂationinGoodsandServices(2019-Present)GoodsTightness(Norm.)ServicesTightness(Norm.)−20246810Inﬂation(YoYGoodsInﬂation(YoYServicesInﬂation(YoY