
# Introduction

This chapter introduces PyPSA-ASEAN and explains how to create PyPSA models in the ASEAN region.

- [Installation](installation.md)
- [Run scenarios](run-scenarios.md)
- [Contributors](contributors.md)

##  Introduction

> All citations are provided in academic publications referenced in the [Index](../index.md)

Southeast Asia is currently experiencing a pivotal moment in the shaping of its future energy landscape. The accelerated economic growth of the region, particularly within the ASEAN-5, where the GDP has increased by 4–8% on an annual basis since the late 1990s, has resulted in a substantial surge in electricity consumption [1]. Concurrently, the region is confronted with mounting climate pressures. In the context of a 3°C global warming scenario, the Asian Development Bank projects a decline in GDP per capita of up to 11% by 2100, thereby posing significant risks to long-term development and energy security [2]. This projected downturn is rooted in cascading climate impacts across critical sectors such as tourism and agriculture [2]. The region's geographic exposure further amplifies these risks. Southeast Asia is among the world's most climate-vulnerable regions, exhibiting a high degree of susceptibility to sea-level rise, intensified tropical storms, and large-scale flooding. These phenomena pose a threat to infrastructure, power systems, and economic resilience [3].

Despite the considerable renewable energy potential in Southeast Asia, which encompasses solar, wind, hydro, geothermal, and biomass sources, the region has witnessed a persistent decline in renewable energy integration. A thorough analysis reveals that the proportion of renewables in total final energy consumption diminished from 38.5% in 2000 to 30.3% in 2015 [4], subsequently declining further to a mere 15.6% of the total primary energy supply in 2022 [5]. This indicates a persistent decline in the significance of renewables within the energy mix. This phenomenon is indicative of a structural lock-in, whereby the prevailing demand for energy continues to be met by fossil fuels rather than clean energy alternatives. Since 2010, the supply of energy in Southeast Asia has been dominated by fossil fuels, which have accounted for approximately 80% of the region's increasing energy demand [6]. This trend is primarily driven by coal, which continues to be a significant energy source due to its abundance, cost-effectiveness, and political support [7]. 

In this context, the integration of regional power grids under the ASEAN Power Grid (APG) has emerged as a key mechanism for reversing the decline in renewable energy. The APG envisions an interconnected regional electricity market designed to facilitate renewable energy deployment by addressing a significant historical impediment: the geographical mismatch between renewable resource abundance and demand centres [8]. The APG facilitates the interconnection of national grids through cross-border transmission, enabling countries with surplus hydropower, solar, or wind energy, such as Laos, to export clean electricity to highly industrialised, fossil-fuel-reliant neighbours [9]. Modelling studies indicate that such interconnection increases system flexibility, enhances variable renewable utilisation, and reduces the need for each country to expand fossil generation capacity to meet rising demand. Modelling studies indicate that such interconnection increases system flexibility, enhances renewable utilisation, and reduces fossil development needs. However, these insights remain largely theoretical in nature unless they are supported by robust modelling tools that can quantify impacts, measure trade-offs, and evaluate long-term outcomes.

In the absence of such a framework, Southeast Asia is susceptible to perpetuating its reliance on fossil fuel infrastructure during a period when rapid decarbonisation is both economically and climatically imperative. It is vital for policymakers to possess the necessary tools that facilitate the exploration of alternative pathways, the comparison of technology strategies, and the evaluation of the implications of regional power integration. Robust, transparent energy system modelling is therefore essential. This is for two main reasons: first, to quantify transition outcomes; and secondly, to identify cost-optimal trajectories. It can also be used to assess cross-border synergies and support evidence-based decision-making. In this context, the development of a high-resolution regional model is of critical importance for the strategic direction of Southeast Asia towards a resilient and low-carbon energy future.


## Problem Statement

Existing modelling efforts are characterised by fragmentation, with efforts limited to individual countries, founded on non-transparent assumptions, exhibiting a lack of high-resolution spatial data or transparency, and often inaccessible to the broader research and policy community [10]. This hinders the capacity to evaluate regional integration strategies, compare transition pathways, or assess policy interventions with methodological transparency.

Open-source energy system models have been posited as a potentially efficacious pathway for regional transition planning, by virtue of the fact that they enable transparent analysis, cross-border collaboration, and evidence-based system design. However, the adoption of these technologies remains limited across Southeast Asia due to the presence of several structural barriers. 

1.  Data availability, consistency, and format standardisation vary widely between countries, making system-level comparison and integration challenging [9]. 
2. Institutional fragmentation and political constraints impede regional cooperation, with progress towards an integrated power network being hindered by non-binding ASEAN agreements and the absence of unified technical and regulatory frameworks. 
3. The capacity for open modelling of the entire region is still in its infancy, with merely two noteworthy initiatives having been identified to date: TZ-APG [11] and the preliminary PyPSA-SEA [12] work. 

These initiatives have exposed a substantial deficit in local modelling infrastructure and expertise.

This gap motivates the research presented in this paper, leading to two central questions:

1. How does a decarbonised transition pathway compare to the baseline in terms of system cost and outcomes?
2. How does regional power integration through the ASEAN Power Grid influence the transition?

To address these questions, this work introduces a fully open-source, high-resolution energy system model of the entire ASEAN region, advancing three key novelties:

1. Full ASEAN-wide modelling, including Timor-Leste
2. High spatial, temporal, and sectoral details
3. Built entirely using open-source tools for transparency and accessibility

A modelling framework that is transparent, high-resolution, and openly accessible is important to represent sectoral interactions, geographic variability, and the system-level impacts of regional interconnection. The establishment of such a platform would facilitate the evaluation of transition pathways on an equal analytical footing by policymakers, as well as the comparison of national and regional strategies, and the identification of least-cost decarbonisation options. The development of PyPSA-ASEAN directly responds to this need by providing the methodological foundation for coordinated, data-driven energy planning, thus supporting a more integrated, secure, and equitable energy future for the region.


