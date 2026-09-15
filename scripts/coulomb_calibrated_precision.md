# Reverse alpha: 80-significant-digit arithmetic archive

These are computed digits from finite-precision measured inputs, not measured digits. See the [input provenance](./coulomb_calibrated_inputs.json), [reproducible calculation](./coulomb_calibrated.py), and [Coulomb exhibit](../coulomb-from-chirality-directed-depletion.html#research-status) for assumptions and scope.

## pi

Formula: `Chudnovsky series; independently checked with Machin identity`

Unit: 1. mathematical constant, not a fitted parameter.

```text
3.1415926535897932384626433832795028841971693993751058209749445923078164062862090E+0
```

## hbar

Formula: `h/(2*pi)`

Unit: J s. fixed SI action unit.

```text
1.0545718176461563912624280033022807447228263300204131224219234705984359127347391E-34
```

## inverse_alpha

Formula: `1/alpha`

Unit: 1. algebraic consequence of supplied inputs.

```text
1.3703599917759012640282278707222856388747882770293902873693170944055273237416711E+2
```

## sqrt_alpha

Formula: `sqrt(alpha)`

Unit: 1. algebraic consequence of supplied inputs.

```text
8.5424543102670440914496789069289043425338720457907491454713775963941047615678975E-2
```

## inverse_sqrt_alpha

Formula: `1/sqrt(alpha)`

Unit: 1. algebraic consequence of supplied inputs.

```text
1.1706237618363558845995168147872770503702483478746271852884795660268390030551696E+1
```

## alpha_squared

Formula: `alpha^2`

Unit: 1. algebraic consequence of supplied inputs.

```text
5.3251354447695785634490000000000000000000000000000000000000000000000000000000000E-5
```

## binding_fraction

Formula: `alpha^2/2`

Unit: 1. quadratic circular model: binding/(mu*c^2).

```text
2.6625677223847892817245000000000000000000000000000000000000000000000000000000000E-5
```

## orbital_to_Compton_period

Formula: `1/alpha^2`

Unit: 1. n=1 circular period divided by h/(mu*c^2); not an event count.

```text
1.8778865070600481799832447088584828003234650799033071563845272362175421769235828E+4
```

## legacy_N

Formula: `2/alpha^2`

Unit: 1. legacy ratio, not a winding integer.

```text
3.7557730141200963599664894177169656006469301598066143127690544724350843538471656E+4
```

## network_kappa

Formula: `4*pi*alpha`

Unit: 1. required weight in the stated network cost; not independently selected.

```text
9.1701236826638077156008688330924432462501840253325677121101474586911448429385112E-2
```

## inverse_network_kappa

Formula: `1/(4*pi*alpha)`

Unit: 1. normalization comparison only.

```text
1.0904978325325186300926589819619849200881212638907912382203669872567129792763447E+1
```

## K_Q

Formula: `alpha*hbar*c`

Unit: J m. algebraic consequence of supplied inputs.

```text
2.3070775507679175960689258194295483335407817925764141330161404878103195769130372E-28
```

## pair_energy_at_1_angstrom

Formula: `K_Q/(1e-10 m)`

Unit: eV. magnitude in the leading 1/r pair law; reference separation only.

```text
1.4399645468602668412582315936019001594931397511420792153533182083657322747555668E+1
```

## pair_force_at_1_angstrom

Formula: `K_Q/(1e-10 m)^2`

Unit: N. magnitude in the leading pair law; reference separation only.

```text
2.3070775507679175960689258194295483335407817925764141330161404878103195769130372E-8
```

## relative_circular_rate_v

Formula: `alpha*c`

Unit: m s^-1. n=1 quadratic relative-motion chart, not motion of light.

```text
2.1876912621441000494000000000000000000000000000000000000000000000000000000000000E+6
```

## connected_phase_per_r_over_c

Formula: `(D/hbar)*(r/c)=alpha`

Unit: rad. reference exposure, not a derived completion duration or propagation time.

```text
7.2973525643000000000000000000000000000000000000000000000000000000000000000000000E-3
```

## connected_phase_degrees

Formula: `alpha*180/pi`

Unit: degree. phase coordinate, distinct from asin(alpha) cone angle.

```text
4.1810750355335867964729843940621875453506735789123619786765425405216192992026522E-1
```

## full_connected_phase_reference_intervals

Formula: `2*pi/alpha`

Unit: 1. stationary phase-rate comparison only; not an integer count of exchanges.

```text
8.6102257658730817818672880461130850289707151497320613594762410657610452006106288E+2
```

## reference_readout_dark_probability

Formula: `sin^2(alpha/2)`

Unit: 1. ideal closed-history readout for c*(I1-I0)=1, st=+1, analysis phase=0; chosen exposure.

```text
1.3312779534804842004314217605728590659954219077531541114407535627804930949007462E-5
```

## cone_angle

Formula: `asin(alpha)`

Unit: rad. legacy cone identification, not a microscopic selection rule.

```text
7.2974173315033574801953005596924823743756816303511844449356030995398578810498480E-3
```

## cone_angle_degrees

Formula: `asin(alpha)*180/pi`

Unit: degree. legacy cone identification.

```text
4.1811121444076192879758061260529070543546140393351516682991108929455608172111084E-1
```

## cone_dilation_minus_one

Formula: `1/sqrt(1-alpha^2)-1`

Unit: 1. algebraic consequence of supplied inputs.

```text
2.6626740661070667396287360306628522387933160764470013467670774094957420152886439E-5
```

## cone_cap_solid_angle

Formula: `2*pi*(1-sqrt(1-alpha^2))`

Unit: sr. algebraic consequence of supplied inputs.

```text
1.6729629114476261625297697132944652632769373755994484598960332745900583040040805E-4
```

## cone_cap_relative_leading_error

Formula: `cap/(pi*alpha^2)-1`

Unit: 1. nonzero correction to the old leading-order equality.

```text
1.3313193087065483672425699181928676629249824943802839558497072814356623927444318E-5
```

## reduced_mass

Formula: `m_light*m_heavy/(m_light+m_heavy)`

Unit: kg. mass inputs; both ends participate.

```text
9.1044252889167705893943168481898187020964735397829715813989242122550870607960392E-31
```

## light_limit_Compton_length

Formula: `hbar/(mass*c)`

Unit: m. mass-calibrated reference length.

```text
3.8615926743523754766874595593809359846441429261617664191685447285519417217553463E-13
```

## light_limit_radius

Formula: `hbar/(mass*c*alpha)`

Unit: m. n=1 quadratic circular model; light_limit uses mass=m_light, two_end uses mass=mu.

```text
5.2917721054674018262542008441642699275770045257394751297319399664132543106806490E-11
```

## light_limit_binding_J

Formula: `mass*c^2*alpha^2/2`

Unit: J. gross binding; no fine, radiative or finite-structure corrections.

```text
2.1798723610794481970184992092421431712253526026502000000000000000000000000000000E-18
```

## light_limit_binding_eV

Formula: `mass*c^2*alpha^2/(2 eV_to_J)`

Unit: eV. gross binding only.

```text
1.3605693122843584030314220705593832616244187484825097006126978631246222506076069E+1
```

## light_limit_binding_mass_equivalent

Formula: `mass*alpha^2/2`

Unit: kg. mass-equivalent of the gross spectral limit; not a new constituent mass.

```text
2.4254351047447815952024439206620550000000000000000000000000000000000000000000000E-35
```

## light_limit_Rydberg

Formula: `mass*c*alpha^2/(2*h)`

Unit: m^-1. gross spectral coefficient; not an independent calibration input.

```text
1.0973731568038493179836040078483540591268717551986677955711048425890872887906265E+7
```

## light_limit_limit_frequency

Formula: `binding/h`

Unit: Hz. gross spectral limit.

```text
3.2898419602144541089992824921150935463992221734178289675970303453548556228309777E+15
```

## light_limit_orbital_period

Formula: `2*pi*hbar/(mass*c^2*alpha^2)`

Unit: s. model recurrence, not an independently observed elementary clock.

```text
1.5198298460738418706683762958987836212689144242487924339857557585473180071071894E-16
```

## light_limit_gross_3_to_2_vacuum_nm

Formula: `1/[R*(1/4-1/9)]`

Unit: nm. uncorrected gross vacuum line; not an air wavelength or resolved observed component.

```text
6.5611227642658373363182651096766159188370246530156841373560513641327916350798910E+2
```

## two_end_Compton_length

Formula: `hbar/(mass*c)`

Unit: m. mass-calibrated reference length.

```text
3.8636957634528894334767221360826371230863222609123794846480496525397184559606271E-13
```

## two_end_radius

Formula: `hbar/(mass*c*alpha)`

Unit: m. n=1 quadratic circular model; light_limit uses mass=m_light, two_end uses mass=mu.

```text
5.2946540946298861196667618655195268801880742657053526688790653756240062278727347E-11
```

## two_end_binding_J

Formula: `mass*c^2*alpha^2/2`

Unit: J. gross binding; no fine, radiative or finite-structure corrections.

```text
2.1786858117019162590149519338368548159833227645559654183474990749851789419598236E-18
```

## two_end_binding_eV

Formula: `mass*c^2*alpha^2/(2 eV_to_J)`

Unit: eV. gross binding only.

```text
1.3598287264136422669954828050730733700010524336207273750239320210839992452166941E+1
```

## two_end_binding_mass_equivalent

Formula: `mass*alpha^2/2`

Unit: kg. mass-equivalent of the gross spectral limit; not a new constituent mass.

```text
2.4241148905133602983236881115986171102751333223466801913553580823497359570906032E-35
```

## two_end_Rydberg

Formula: `mass*c*alpha^2/(2*h)`

Unit: m^-1. gross spectral coefficient; not an independent calibration input.

```text
1.0967758340158852161776037620149972799340786861183737734897556902227993439326544E+7
```

## two_end_limit_frequency

Formula: `binding/h`

Unit: Hz. gross spectral limit.

```text
3.2880512315462224000374519636452306741475152727683775251722909619137958295835784E+15
```

## two_end_orbital_period

Formula: `2*pi*hbar/(mass*c^2*alpha^2)`

Unit: s. model recurrence, not an independently observed elementary clock.

```text
1.5206575712777824330017341511098898347891573856858952395936288426632457461433229E-16
```

## two_end_gross_3_to_2_vacuum_nm

Formula: `1/[R*(1/4-1/9)]`

Unit: nm. uncorrected gross vacuum line; not an air wavelength or resolved observed component.

```text
6.5646960634033429878548668717023931442014378234122511172391299310966198023186513E+2
```

## equal_end_exact_radius

Formula: `2*hbar*sqrt(1-alpha^2/4)/(m_light*c*alpha)`

Unit: m. full two-end kinetic chart, equal masses and n=1, inherited 1/r potential.

```text
1.0583473762192323518876112452320243723839301582710974432190989696514749354673216E-10
```

## equal_end_exact_binding_eV

Formula: `2*m_light*c^2*[1-sqrt(1-alpha^2/4)]/eV_to_J`

Unit: eV. conditional model result, not a fine-structure calculation.

```text
6.8028692028720967254012736611740651143890294682668451344654082815884505984656653E+0
```

## equal_end_binding_fractional_correction

Formula: `binding_exact/binding_quadratic-1`

Unit: 1. nonlinear kinetic correction within this trial model.

```text
3.3282318071243092619962186026877375029333098082016782847031539192449595756300539E-6
```

## gravitational_pair_dimensionless

Formula: `G*m_light*m_heavy/(hbar*c)`

Unit: 1. G and both masses supplied independently of alpha.

```text
3.2165895109146318363348133709287179648931927714924626624605890511626894009775911E-42
```

## source_Ks_opposite

Formula: `mu*c^2/2 * [sqrt(alpha^4-gamma_g^4*(1-u^2))-gamma_g^2*u]`

Unit: J. conditional reverse circular map C_*=K_Q, n=1; u remains unspecified.

```text
2.1786858117019162590149519338368548159833227645559654183474990749851789419598240E-18
```

## source_Ks_relative_shift_opposite

Formula: `(required_Ks - mu*c^2*alpha^2/2)/(mu*c^2*alpha^2/2)`

Unit: 1. stable small-difference formula; far below physical input precision.

```text
1.9429455248670631213118283129845093777151830630046197100689541670544358578338689E-79
```

## source_Ks_absolute_shift_opposite

Formula: `base_energy * relative_shift`

Unit: J. mathematical correction, not experimentally resolved.

```text
4.2330678479376031379279832633319048272997487879078450698603432833022701176118677E-97
```

## source_Ks_perpendicular

Formula: `mu*c^2/2 * [sqrt(alpha^4-gamma_g^4*(1-u^2))-gamma_g^2*u]`

Unit: J. conditional reverse circular map C_*=K_Q, n=1; u remains unspecified.

```text
2.1786858117019162590149519338368548159833227645559654183474990749851789419598236E-18
```

## source_Ks_relative_shift_perpendicular

Formula: `(required_Ks - mu*c^2*alpha^2/2)/(mu*c^2*alpha^2/2)`

Unit: 1. stable small-difference formula; far below physical input precision.

```text
-1.8875186563004736989541850571392287275888722531903129025815027421823556485636952E-158
```

## source_Ks_absolute_shift_perpendicular

Formula: `base_energy * relative_shift`

Unit: J. mathematical correction, not experimentally resolved.

```text
-4.1123101158045078345415258558201094563807756001908055661211236390864574460170732E-176
```

## source_Ks_aligned

Formula: `mu*c^2/2 * [sqrt(alpha^4-gamma_g^4*(1-u^2))-gamma_g^2*u]`

Unit: J. conditional reverse circular map C_*=K_Q, n=1; u remains unspecified.

```text
2.1786858117019162590149519338368548159833227645559654183474990749851789419598232E-18
```

## source_Ks_relative_shift_aligned

Formula: `(required_Ks - mu*c^2*alpha^2/2)/(mu*c^2*alpha^2/2)`

Unit: 1. stable small-difference formula; far below physical input precision.

```text
-1.9429455248670631213118283129845093777151830630046197100689541670544358578338689E-79
```

## source_Ks_absolute_shift_aligned

Formula: `base_energy * relative_shift`

Unit: J. mathematical correction, not experimentally resolved.

```text
-4.2330678479376031379279832633319048272997487879078450698603432833022701176118677E-97
```

## Planck_mass

Formula: `sqrt(hbar*c/G)`

Unit: kg. uses measured G; not derived from alpha.

```text
2.1764343427178982139279149190241470411173913749165146773123128134687551077270348E-8
```

## Planck_frequency

Formula: `M_Pl*c^2/h`

Unit: Hz. uses measured G; cycle frequency, not angular frequency.

```text
2.9520991966835316844042440449607740132268257526791444545297815867877759556249623E+42
```

## light_mass_over_Planck

Formula: `m_light/M_Pl`

Unit: 1. uses measured mass and G, not a result of alpha.

```text
4.1854622191470934565925602304292742773305344287303868092806003273252294077298636E-23
```

## alpha_G_light

Formula: `G*m_light^2/(hbar*c)`

Unit: 1. uses mass and G.

```text
1.7518093987907712171482191126859670679426084313794860411098105994705800591259680E-45
```

## rate_zeta_0p25

Formula: `sqrt(alpha/zeta)*f_Pl`

Unit: Hz. conditional normalization example; not an independently established completion clock.

```text
5.0436345014090227182616549994817550194407066457109694191444661466902034223902972E+41
```

## rate_zeta_1

Formula: `sqrt(alpha/zeta)*f_Pl`

Unit: Hz. conditional normalization example; not an independently established completion clock.

```text
2.5218172507045113591308274997408775097203533228554847095722330733451017111951486E+41
```

## rate_zeta_4

Formula: `sqrt(alpha/zeta)*f_Pl`

Unit: Hz. conditional normalization example; not an independently established completion clock.

```text
1.2609086253522556795654137498704387548601766614277423547861165366725508555975743E+41
```

## reference_rate_period_in_Planck_times

Formula: `2*pi/sqrt(alpha)`

Unit: 1. requires zeta=1; no count or actual period established.

```text
7.3552460206054868297318887474146507826671948027585648608151582566339148309810850E+1
```

## spectral_magnitude_product_limit

Formula: `8*alpha/pi`

Unit: 1. magnitude only of (Omega*a/c)*(r/ell); actual required product is negative in this trial.

```text
1.8582555713482607984324375084720833534891882572943831016340189068984974663122899E-2
```

## spectral_required_product_N3

Formula: `-2*alpha/[N*tan(pi/(4*N))]`

Unit: 1. no solution with positive frequency, lengths and c under the original sign assignment.

```text
-1.8156060353807088377863036235079020826945726980219128948068058390006267206354178E-2
```

## spectral_required_product_N10

Formula: `-2*alpha/[N*tan(pi/(4*N))]`

Unit: 1. no solution with positive frequency, lengths and c under the original sign assignment.

```text
-1.8544331142809057184288286037399455345391472199143665884235192338635342618417913E-2
```

## spectral_required_product_N100

Formula: `-2*alpha/[N*tan(pi/(4*N))]`

Unit: 1. no solution with positive frequency, lengths and c under the original sign assignment.

```text
-1.8582173623424545173608930655670337372624102264415942874230373911616595274025491E-2
```

## relative_input_alpha_uncertainty

Formula: `u(alpha)/alpha`

Unit: 1. measurement precision, irrespective of arithmetic precision.

```text
1.5073959909534913904310506577945142027622671047323293161062488038460800561158382E-10
```

## inverse_alpha_first_order_uncertainty

Formula: `u(alpha)/alpha^2`

Unit: 1. first-order standard uncertainty from alpha alone.

```text
2.0656751577660529979815691797443310803558115878936378720229799598392963946159411E-8
```

## legacy_N_first_order_uncertainty

Formula: `4*u(alpha)/alpha^3`

Unit: 1. first-order standard uncertainty from alpha alone.

```text
1.1322874368831887729607742832211460129277786320886082239818360202998456595297662E-5
```

## inverse_alpha_printed_reference_difference

Formula: `1/(printed alpha) - separately printed inverse alpha`

Unit: 1. small independent-rounding discrepancy; not physical evidence.

```text
5.9012640282278707222856388747882770293902873693170944055273237416711174923435787E-10
```

## Rydberg_printed_reference_relative_difference

Formula: `R_from_printed_mass_and_alpha/R_CODATA-1`

Unit: 1. correlated adjustment and finite printed rounding; not an independent precision prediction.

```text
-1.0799136048473867422455709511296923462348029188220702203946644036200754340756112E-11
```

## Rydberg_printed_rounding_allowance

Formula: `half-last-digit relative allowances for mass + 2*alpha + R`

Unit: 1. conservative first-order input/output rounding allowance, not statistical uncertainty.

```text
1.9238008791605360168764550652644425765596578182521634826579165636354669123781647E-11
```
