import os
import numpy as np
import astropy.units as u

def ensure_dirs():
    """Ensure Figures and Report directories exist."""
    for folder in ['Figures', 'Report']:
        if not os.path.exists(folder):
            os.makedirs(folder)
            print(f"'{folder}' directory was created!")

def generate_report_name(event, report_name=None):
    """Generate a default report name if none is provided."""
    body_name = event.body.shortname.replace('(','').replace(')','').replace('/','').replace(' ','_')
    if report_name is None:
        report_name = f'Report_{event.tca.iso[:10].replace("-","")}_{body_name}.tex'
    return report_name

def write_header(f, event, author=None):
    """Write LaTeX header."""
    body_name = event.body.shortname.replace('_','\_')
    header = r'''\documentclass{article} '''
    header += '\n\\usepackage[utf8]{inputenc}'
    header += '\n\\usepackage[left=2cm, right=2cm, top=2.5cm]{geometry}'
    header += '\n\\usepackage{helvet}'
    header += '\n\\usepackage{tabularx, graphicx, xcolor}'
    header += '\n\\usepackage{natbib}'
    header += '\n\\bibpunct{(}{)}{;}{a}{}{,} % to follow the A&A style'
    header += '\n\\begin{document} \n'
    header += f'\n\\title{{\\textbf{{Stellar occultation by {body_name}\\\on {event.tca.datetime.strftime("%d %B %Y")} }}}}'
    header += f'\n\\author{{\\textbf{{{author if author else "Author name"}}}}}'
    header += '\n\\date{\\today}'
    header += '\n\\maketitle \n'
    header += '\\begin{center} \nPreliminary report generated automatically with SORA v0.3.2-dev \n\\end{center} \n'
    header += '\n\\clearpage\n\n'
    f.write(header)

def write_prediction(f, event):
    """Write prediction section."""
    f.write(r'\section{Prediction}' + '\n')
    body_name = event.body.shortname.replace('_','\_')
    map_name = 'Figures/{}_{}'.format(event.tca.iso[:10].replace('-',''), body_name.replace('(','').replace(')','').replace('/','').replace(' ','_').replace('-',''))
    pred_map = event.plot_occ_map(nameimg=map_name, site_name=False, sscale=1);

    prediction = '\n\\begin{figure}[h]'
    prediction += '\n\t\centering'
    prediction += f'\n\t\includegraphics[width=0.9\hsize]{{{map_name}.png}}'
    prediction += '\n\t\caption{Prediction map for the occultation event analyzed in this report. The blue lines show the limits corresponding to the occulting body shadow projected on the Earth. The red dashed lines extend these limits considering the 1-sigma error bar.}'    
    prediction += '\n\t\label{fig:pred_map}'
    prediction += '\n\end{figure}\n\n'
    prediction += '\clearpage\n\n'
    f.write(prediction)

def write_occultation_params(f, event):
    """Write occultation parameters table."""
    occ_circums = r'\begin{table}[ht]' + '\n\centering\n\small\n\caption{Event circumstances.}\n'
    occ_circums += '\t\\begin{tabular}{c c} \hline \hline\n'
    occ_circums += f'\tEpoch & {event.tca.iso} UTC \\\ \n'
    occ_circums += f'\tStar position (ICRS) & {event.star.ra.to_string()}, {event.star.dec.to_string()} \\\ \n'
    occ_circums += f'\tClosest Approach & {round(event.ca.value,3)} arcsec \\\ \n'
    occ_circums += f'\tPosition Angle & {round(event.pa.value,2)} deg \\\ \n'
    occ_circums += f'\tShadow velocity & {round(event.vel.value,2)} km/s \\\ \n'
    occ_circums += f'\tGeocentric distance & {round(event.dist.value,2)} au \\\ \n'
    occ_circums += '\\hline\\end{tabular}\n\\end{table}\n\\clearpage\n'
    f.write(occ_circums)

def write_star_params(f, event, report_name="report.tex"):
    """
    Writes the stellar parameters in LaTeX format.
    """
    import astropy.units as u

    star = getattr(event, 'star', None)
    tca = getattr(event, 'tca', None)
    dist = getattr(event, 'dist', None)
    ref_center = getattr(event, '_reference_center', None)

    star_params = r'''\begin{table}[ht]'''
    star_params += '\n\centering'
    star_params += '\n\small'
    star_params += '\n\caption{Parameters of the target star.}'
    star_params += '\n\\vspace{2mm}'
    star_params += '\n\t\\begin{tabular}{c c} \hline \hline'

    star_params += '\n\tStar source ID \t & {} \\\ '.format(star.code)
    star_params += '\n \\hline'
    star_params += '\n\tStellar catalogue \t & {} \\\ '.format(star._catalogue)
    star_params += '\n \\hline'
    star_params += '\n\tPosition in catalog (ICRS) \t & RA: {} $\pm$ {:.5f} \\\ '.format(
        star.coord.ra.to_string(u.hourangle, sep='hms', precision=5), 
        star.errors['ra']
    )
    star_params += '\n\t\t & DEC: {} $\pm$ {:.4f} \\\ '.format(
        star.coord.dec.to_string(u.deg, sep='dms', precision=4), 
        star.errors['dec']
    )
    star_params += '\n \\hline'
    star_params += '\n\tProper motion \t & pmRA: {:.3f} $\pm$ {:.3f} mas/yr\\\ '.format(
        star.pmra.value, star.errors['pmra'].value
    )
    star_params += '\n\t & pmDEC: {:.3f} $\pm$ {:.3f} mas/yr\\\ '.format(
        star.pmdec.value, star.errors['pmdec'].value
    )
    star_params += '\n \\hline'
    star_params += '\n\tParallax \t & {:.4f} $\pm$ {:.4f} mas\\\ '.format(
        star.parallax.value, star.errors['parallax'].value
    )
    star_params += '\n \\hline'
    star_params += '\n\tRadial velocity \t & {:.2f} $\pm$ {:.2f} km/s\\\ '.format(
        star.rad_vel.value, star.errors['rad_vel'].value
    )
    star_params += '\n \\hline'
    coord = star.get_position(tca, observer=ref_center)
    try:
        error_star = star.error_at(tca)
    except:
        error_star = [0, 0]*u.mas
    star_params += '\n\tGeocentric star coordinate (epoch) & RA: {} $\pm$ {:.5f} \\\ '.format(
        coord.ra.to_string(u.hourangle, sep='hms', precision=5), error_star[0]
    )
    star_params += '\n\t & DEC: {} $\pm$ {:.4f} \\\ '.format(
        coord.dec.to_string(u.deg, sep='dms', precision=4), error_star[1]
    )
    star_params += '\n \\hline'

    # Magnitudes
    mag_out = ['& {}: {:6.3f} \\\ '.format(mag, star.mag[mag]) for mag in star.mag]
    out_mag = []
    for i, mag in enumerate(mag_out):
        if i % len(mag) == 0:
            out_mag.append([])
        out_mag[-1].append(mag)
    star_params += '\n\tMagnitudes'
    for mags in out_mag[0]:
        star_params += mags
    star_params += '\n \\hline'
    
    # Diâmetros Angulares (Kervella e van Belle)
    kerv = event.star.kervella()
    vanb = event.star.van_belle()
        
    if kerv:
        star_params += '\n\tAngular diameter from Kervella et al. (2004) &   ' + ','.join([' {}: {:.4f}'.format(k, v) for k, v in kerv.items()])   
        star_params += ' \\\ '

    else:
        star_params += "\n\tAngular diameter (Kervella et al. (2004) & Mag B and V not furnished \\\ \hline"
    star_params += '\n \\hline'

    if vanb:
        star_params += '\n\tAngular diameter from van Belle (1999)'
        for key, value in vanb.items():
            star_params += ' &  {}:'.format(key)
            star_params += ','.join([' {}: {:.4f} '.format(k, v) for k, v in value.items()])
            star_params += ' \\\ '
    else:
        star_params += "\n\tAngular diameter (van Belle, 1999) & Mag B and V not furnished \\\ \hline"
    star_params += '\n \\hline'
    # Diâmetros Aparentes
    if kerv:
        star_params += '\n\tApparent diameter from Kervella et al. (2004)'
        for k, v in kerv.items():
            star_diam = event.star.apparent_diameter(distance=event.dist, band=k, mode='kervella', verbose=False);
            star_params += '& {}: {:.3f} \\\ \n '.format(k, star_diam)
    else:
        star_params += "\n\tApparent diameter (Kervella+ 2004) & Could not be calculated (missing magnitudes) \\\\"
    star_params += '\n \\hline'

    if vanb:
        star_params += '\n\tApparent diameter from van Belle (1999)'
        for k, v in vanb.items():
            for h, j in v.items():
                star_diam = event.star.apparent_diameter(distance=event.dist, band=h, mode='van_belle', verbose=False, star_type=k);
                star_params += '& {}: {}: {:.3f} \\\ \n '.format(k, h, star_diam)
    else:
        star_params += "\n\tApparent diameter (van Belle 1999) & Could not be calculated (missing magnitudes) \\\\"

        
    star_params += '\hline\n\t\end{tabular}'
    star_params += '\n\label{tab:star_params}'
    star_params += '\n\end{table}\n\n\clearpage \n'

    f.write(star_params)


def write_observational_circumstances(f, event):
    """Gera a tabela de circunstâncias observacionais no formato de múltiplas linhas."""
    header = r"""\begin{table}[ht]
            \centering
            \caption{Observational circumstances of the sites involved in this campaign.}
            \label{tab:obs_circumstances}
            \small
            \begin{tabular}{l l l l l l} \hline \hline
             & Latitude & Telescope (mm) & Exposure (s) & & \\
            Site & Longitude & Camera & Cycle (s) & Status & Observers \\
             & Altitude (m) & Filter & & & \\ \hline
            """

    footer = r"""\end{tabular}
            \end{table}
            \clearpage
            """
    
    f.write(header)
    
    for i, chord in enumerate(event.chords):
        observer = event.chords[i].observer
        lc = event.chords[i].lightcurve

        # Extração segura de todos os dados
        site_name = getattr(lc, 'name', 'N/A')
        latitude = observer.lat.to_string(sep=' ', precision=3) if hasattr(observer, 'lat') else 'N/A'
        telescope = getattr(observer, 'telescope', 'N/A')
        exposure = f"{lc.exptime:.4f}" if hasattr(lc, 'exptime') else 'N/A'
        longitude = observer.lon.to_string(sep=' ', precision=3) if hasattr(observer, 'lon') else 'N/A'
        camera = getattr(observer, 'camera', 'N/A')
        cycle = f"{lc.cycle:.4f}" if hasattr(lc, 'cycle') else 'N/A'
        status = event.chords[i].status() if hasattr(event.chords[i], 'status') else 'N/A'
        observers = ", ".join(observer.observers) if hasattr(observer, 'observers') and isinstance(observer.observers, list) else getattr(observer, 'observers', 'N/A')
        altitude_val = observer.height.to('m').value if hasattr(observer, 'height') else 'N/A'
        altitude = f"{altitude_val:.2f}" if isinstance(altitude_val, (int, float)) else altitude_val
        filter_name = getattr(observer, 'filter', 'N/A')

        # Escrita das linhas da tabela para a corda atual
        f.write(f'{site_name} & {latitude} & {telescope} & {exposure} & & \\\\ \n')
        f.write(f' & {longitude} & {camera} & {cycle} & {status} & {observers} \\\\ \n')
        f.write(f' & {altitude} & {filter_name} & & & \\\\ \hline \n')
        
    f.write(footer)

def generate_scientific_text(lc, event, chord_index):
    """
    Gera uma descrição textual cientificamente aprimorada de uma curva de luz de ocultação.

    """
    
    text_parts = []
    
    if hasattr(lc, "initial_time") and hasattr(lc, "end_time"):
        duration_min = (lc.end_time - lc.initial_time).to("min").value
        observer_name = event.chords[chord_index].observer.name
        obs_date = event.tca.iso[:10]
        start_time_utc = lc.initial_time.iso[11:]
        end_time_utc = lc.end_time.iso[11:]        
        flux = getattr(lc, "flux", None)

        if flux is not None and len(flux) > 0:
            n_points = len(flux)
        else:
            n_points = 'NaN'


        p1 = (
            f"The stellar occultation was observed from the {observer_name} on {obs_date}. "
            f"Data acquisition started at {start_time_utc} UTC and ended at {end_time_utc} UTC, "
            f"spanning a total of {duration_min:.2f} minutes. "
        )
        
        cycle = getattr(lc, "cycle", getattr(lc, "exptime", None))
        if lc.exptime and cycle:
            p1 += (
                f"A total of {n_points} data points were collected using an exposure time of {lc.exptime:.4f} s "
                f"and a cycle time of {cycle:.4f} s. "
            )
            
        if hasattr(lc, "central_lambda") and hasattr(lc, "delta_lambda"):
            p1 += (
                f"A photometric filter centered at {lc.central_lambda:.2f} μm "
                f"(FWHM = {lc.delta_lambda:.2f} μm) was used. "
            )
        
        text_parts.append(p1)

    p2_parts = []
    if hasattr(lc, "dist") and hasattr(lc, "vel"):
        p2_parts.append(
            f"At the time of the event, the object was at a geocentric distance of {lc.dist:.4f} au, "
            f"with a shadow velocity of {lc.vel:.3f} km s$^{{-1}}$. "
        )

        fresnel_t = lc.fresnel_scale / lc.vel if hasattr(lc, "fresnel_scale") else 0
        dstar_t = lc.d_star / lc.vel if hasattr(lc, "d_star") else 0
        inst_response_km = lc.exptime * lc.vel if hasattr(lc, "exptime") else 0
        
        resolution_text = []
        if fresnel_t > 0:
            resolution_text.append(f"a Fresnel scale of {lc.fresnel_scale:.3f} km ({fresnel_t:.3f} s)")
        if dstar_t > 0:
            resolution_text.append(f"a projected stellar diameter of {lc.d_star:.3f} km ({dstar_t:.3f} s)")
        if inst_response_km > 0:
            resolution_text.append(f"a instrumental response of {inst_response_km:.3f} km due to the exposure time")
        
        if resolution_text:
            p2_parts.append(
                "The effective spatial resolution is mainly determined by " + ", ".join(resolution_text) + ". "
            )

        if hasattr(lc, "model_resolution"):
             model_res_km = lc.model_resolution * lc.vel
             p2_parts.append(
                 f"The final model resolution was {lc.model_resolution:.3f} s "
                 f"({model_res_km:.3f} km). "
            )

    if p2_parts:
        text_parts.append("".join(p2_parts))

    p3_parts = []
    if hasattr(lc, "baseflux") and hasattr(lc, "bottomflux"):
        p3_parts.append(
            f"The light curve was normalized to a baseline flux of {lc.baseflux:.3f}. "
            f"During the event, the flux dropped to a minimum of {lc.bottomflux:.3f}. "
        )

    if hasattr(lc, "lc_sigma") and hasattr(lc, 'baseflux'):
        mean_sigma = lc.lc_sigma.mean()
        snr = lc.baseflux / mean_sigma if mean_sigma > 0 else float('inf')
        p3_parts.append(
            f"The mean photometric uncertainty of the out-of-event data points was {mean_sigma:.3f}, "
            f"which corresponds to an average signal-to-noise ratio (SNR) of approximately {snr:.1f} per point. "
        )
    
    if p3_parts:
        text_parts.append("".join(p3_parts))
        
    return "\n\n".join(text_parts)


def write_light_curves(f, event):
    """Write detailed light curve section."""
    f.write(r'\section{Light curve modeling}' + '\n')
    for i, chord in enumerate(event.chords):
        f.write(f'\n\\subsection{{{event.chords[i].observer.name}}}\n')
        lc = event.chords[i].lightcurve
        
                    
        scientific_text = generate_scientific_text(lc, event, i)
        f.write(scientific_text)
        if hasattr(lc, 'flux') and lc.flux is not None:    
            lc_name = f'Figures/LC_full_{lc.name}.pdf'
            model_name = f'Figures/LC_modeled_{lc.name}.pdf'
            
            _report_lc_plot(event, chord=event.chords[i], filename=lc_name)
            _report_model_plot(event, chord=event.chords[i], filename=model_name)
            
            f.write(f'''
            
            \\begin{{figure}}[h]
                \\centering
                \\includegraphics[width=\\hsize]{{{lc_name}}}
                \\caption{{Full observed light curve.}}
            \\end{{figure}}

            \\begin{{figure}}[h]
                \\centering
                \\includegraphics[width=\\hsize]{{{model_name}}}
                \\caption{{Modeled light curve.}}
            \\end{{figure}}
            ''')
        f.write('\n\\clearpage\n')

def write_ellipse_fit(f, event):
    """Write ellipse fit section."""
    if hasattr(event, 'fitted_params'):
        ellipse_fit = '\n \clearpage \n '
        ellipse_fit += '\n\section{Bi-dimensional modeling}\n '

        new_astrometry = event.new_astrometric_position(verbose=False)
        apparent_b = event.fitted_params['equatorial_radius'][0]*(1.0-event.fitted_params['oblateness'][0])
        equiv_r = np.sqrt(event.fitted_params['equatorial_radius'][0]**2*(1.0-event.fitted_params['oblateness'][0]))

        if event.fitted_params['oblateness'][0] == 0.000000:
            ellipse_fit += 'A circle with radius {:.2f} km was fitted to the positive chord(s). '.format(
                event.fitted_params['equatorial_radius'][0]
            )
            ellipse_fit += 'The center solution for the circle is f$_c = {:.2f} \pm {:.2f}$~km and g$_c = {:.2f} \pm {:.2f}$~km, with uncertainties at $1\sigma$ level. '.format(
                event.fitted_params['center_f'][0],
                event.fitted_params['center_f'][1],
                event.fitted_params['center_g'][0],
                event.fitted_params['center_g'][1]
            )
            ellipse_fit += 'These values correspond to an {}. '.format(
                new_astrometry.split('\n')[0].replace('+/-', '$ \pm $')
            )

            ellipse_fit += 'In milli-arcseconds, the offset is {}. '.format(
                new_astrometry.split('\n')[1].replace('+/-', '$ \pm $').replace('Ephemeris offset (mas): ', '').replace('da_cos_dec ', 'da\_cos\_dec').replace('d_dec ', 'd\_dec')
            )
            ellipse_fit += '{} {}. \n '.format(
                new_astrometry.split('\n')[3].replace('+/-', '$ \pm $'), 
                new_astrometry.split('\n')[4].replace('+/-', '$ \pm $')
            )

        else:
            ellipse_fit += "The projected ellipse fitted to the extremities of the positive chords has a semi-major axis of $a' = {:.2f} \pm {:.2f}$~km. ".format(
                event.fitted_params['equatorial_radius'][0],
                event.fitted_params['equatorial_radius'][1]
            )
            ellipse_fit += "The apparent oblateness is $\epsilon' = {:.3f} \pm {:.3f}$, resulting in a semi-minor axis of $b' = {:.3f} \pm $~km. ".format(
                event.fitted_params['oblateness'][0],
                event.fitted_params['oblateness'][1],
                apparent_b
            )
            ellipse_fit += "The equivalent radius can be found using the equation R$_{{equiv}} = \sqrt{{{{a'}}^{{2}}\,(1 - \epsilon')}}$, resulting in R$_{{equiv}} = {:.2f}$~km. ".format(
                equiv_r
            )

            ellipse_fit += 'The center solution for the circle is f$_c = {:.2f} \pm {:.2f}$~km and g$_c = {:.2f} \pm {:.2f}$~km, with uncertainties at $1\sigma$ level. '.format(
                event.fitted_params['center_f'][0], 
                event.fitted_params['center_f'][1],
                event.fitted_params['center_g'][0],
                event.fitted_params['center_g'][1]
            )
            ellipse_fit += 'These values correspond to an {}. '.format(
                new_astrometry.split('\n')[0].replace('+/-', '$ \pm $')
            )

            ellipse_fit += 'In milli-arcseconds, the offset is {}. '.format(
                new_astrometry.split('\n')[1].replace('+/-', '$ \pm $').replace('Ephemeris offset (mas): ', '').replace('da_cos_dec ','da\_cos\_dec').replace('d_dec ','d\_dec')
            )
            ellipse_fit += '{} {}. \n '.format(
                new_astrometry.split('\n')[3].replace('+/-', '$ \pm $'), 
                new_astrometry.split('\n')[4].replace('+/-', '$ \pm $')
            )

        if not np.isnan(event.body.H):
            H_sun = -26.74
            geometric_albedo = (10**(0.4*(H_sun - event.body.H.value))) * ((u.au.to('km')/equiv_r)**2)
            phys_par = event.body.meta_sbdb['phys_par']
            ellipse_fit += "Using the area of the projected limb (R$_{{equiv}}$), the body absolute magnitude (H$_{{V}} = {:.3f}$ mag [ref. {}]), and the sun absolute albedo (H$_{{\odot}} = {}$), we determine the instantaneous geometric albedo as $p_V = {:.3f}$ (or {:.1f}\%). ".format(
                event.body.H.value, phys_par['H_ref'],  
                -26.74, 
                (10**(0.4*(H_sun - event.body.H.value))) * ((u.au.to('km')/equiv_r)**2), 
                100*(10**(0.4*(H_sun - event.body.H.value))) * ((u.au.to('km')/equiv_r)**2)
            )
        else:
            ellipse_fit += 'The geometric albedo (V) was not calculated as the absolute magnitude (H) is unknown. \n'

        ellipse_fit += '\\textbf{The $\chi^2$ figures from the ellipse fit were not generated automatically. Please insert by hand. }'  

        ellipse_fit += '\\begin{table}[ht] \n '
        ellipse_fit += '\centering \n '
        ellipse_fit += '\caption{Parameters obtained during the bi-dimensional fit procedure.} \n '
        ellipse_fit += '\\vspace{2mm} \n '
        ellipse_fit += '\\begin{tabular}{c c} \hline \hline  \n '
        ellipse_fit += 'Parameter & Value \\\ \n \hline \n '

        for param, value in event.fitted_params.items():
            ellipse_fit += '{} & ${:.3f} \pm {:.3f}$ \\\ \n '.format(
                param.replace('_', '\_'), value[0], value[1]
            )

        ellipse_fit += '\n\hline\n\end{tabular} \n'
        ellipse_fit += '\end{table} \n'

        ellipse_fit += '\\begin{table}[ht] \n '
        ellipse_fit += '\centering \n '
        ellipse_fit += '\caption{$\chi^2$ results from the bi-dimensional fit procedure.} \n '
        ellipse_fit += '\\vspace{2mm} \n '
        ellipse_fit += '\\begin{tabular}{c c} \hline \hline  \n '
        ellipse_fit += 'Parameter & Value \\\ \n \hline \n '
        ellipse_fit += '$\chi^2_{{min}}$      &  {:.3f} \\\ \n'.format(
            event.chi2_params['chi2_min']
        )
        ellipse_fit += 'N$_{{points}}$         &  {} \\\ \n'.format(
            event.chi2_params['npts']
        )
        ellipse_fit += '$\chi^2_{{min}} \\textit{{pdf}}$   &  {:.3f} \\\ \n'.format(
            event.chi2_params['chi2_min']/(event.chi2_params['npts'] - event.chi2_params['nparam'])
        )
        ellipse_fit += 'Radial dispersion     &  {:.3f} $\pm$ {:.3f}~km \\\ \n'.format(
            event.chi2_params['radial_dispersion'].mean(), 
            event.chi2_params['radial_dispersion'].std(ddof=1)
        )
        ellipse_fit += 'Radial error          &  {:.3f} $\pm$ {:.3f}~km \\\ \n'.format(
            event.chi2_params['radial_error'].mean(), 
            event.chi2_params['radial_error'].std(ddof=1)
        )
        ellipse_fit += '\n\hline\n\end{tabular} \n'
        ellipse_fit += '\end{table} \n\clearpage \n'



        #try:
        _projected_chords_plot(event)
        ellipse_fit += '\\begin{figure}[ht] \n'
        ellipse_fit += '\centering \n'
        ellipse_fit += '\includegraphics[width=\hsize]{Figures/projected_chords_plot.pdf} \n '
        ellipse_fit += '\caption{Occultation chords and best-fitted ellipse plotted onto sky plane.} \n '   
        ellipse_fit += '\end{figure} \n'

    f.write(ellipse_fit)

def finalize_document(f):
    """Close LaTeX document."""
    f.write('\n\\end{document}\n')

def to_report(event, pdf=False, author=None, report_name=None):
    """Generate full LaTeX report for an occultation event."""
    ensure_dirs()
    report_name = generate_report_name(event, report_name)
    report_path = os.path.join('Report', report_name)
    
    with open(report_path, 'w') as f:
        write_header(f, event, author)
        write_prediction(f, event)
        write_occultation_params(f, event)
        write_star_params(f, event, report_name=report_name)
        write_observational_circumstances(f, event)
        write_light_curves(f, event)
        write_ellipse_fit(f, event)
        finalize_document(f)
    
    if pdf:
        os.system(f"pdflatex {report_path}")
        base_name = os.path.splitext(report_name)[0]
        os.system(f"mv {base_name}.* Report/")

    return report_path

def _report_lc_plot(event, chord, filename):        
        print("Generating the observed light curves...")
        import matplotlib.pyplot as pl
       
        ''' 
        
        Generates the graph for the light curves with 
        the pattern chosen by me. It fits well in the pdf 
        report.
        C.L.P.
        
        '''
        
        pl.figure(figsize=(8,3))
        pl.plot(chord.lightcurve.time, chord.lightcurve.flux, 'k.-', label='data', zorder=1)

        pl.axhline(y=1, linewidth=1, color='grey', linestyle=':', zorder=0)
        pl.axhline(y=0, linewidth=1, color='grey', linestyle=':', zorder=0)

        pl.xlabel('Instants from {} UTC (seconds)'.format(chord.lightcurve.tref.iso))
        pl.ylabel('Normalized flux ratio')
        pl.legend(ncols=3)
        pl.title('Occultation by {}  @ {}'.format(event.body.shortname,
                                                  chord.observer.name),
                 fontsize=8)
        pl.ylim(0 - chord.lightcurve.flux.std(ddof=1), np.median(chord.lightcurve.flux)+chord.lightcurve.flux.std(ddof=1)*3)
        pl.legend(loc=1, ncol=3)
        pl.tight_layout()
        pl.savefig(filename, dpi=400)
        pl.clf() 
        print("\tFigure '{}' was generated!".format(
            filename
        ))
    
def _report_model_plot(event, chord, filename):
        print("Generating the modeled light curves...")
        import matplotlib.pyplot as pl
        
        pl.figure(figsize=(8,3))
        pl.plot(chord.lightcurve.time, chord.lightcurve.flux, 'k.-', label='data', zorder=1)
        pl.plot(chord.lightcurve.time, chord.lightcurve.model, 'r-', label='model', zorder=2)
        pl.scatter(chord.lightcurve.time, chord.lightcurve.model, s=20, color='r', zorder=2)
        try:
            pl.plot(chord.lightcurve.time_model, chord.lightcurve.model_geometric, 'c-', label='geometric', zorder=0)
        except:
            pass
        #pl.plot(chord.lightcurve.time_model, chord.lightcurve.model_fresnel, 'c-', label='fresnel', zorder=0)
        #pl.plot(chord.lightcurve.time_model, chord.lightcurve.model_star, 'c-', label='star', zorder=0)

        pl.axhline(y=1, linewidth=1, color='grey', linestyle=':', zorder=0)
        pl.axhline(y=0, linewidth=1, color='grey', linestyle=':', zorder=0)

        pl.xlabel('Instants from {} UTC (seconds)'.format(chord.lightcurve.tref.iso))
        pl.ylabel('Normalized flux ratio')
        pl.legend(ncols=3)
        pl.title('Occultation by {}  @ {}'.format(event.body.shortname, 
                                                  chord.observer.name
                                                 ),
                 fontsize=8)
        if chord.status() == 'positive':
            pl.xlim((chord.lightcurve.immersion - chord.lightcurve.tref).sec - chord.lightcurve.exptime*40,
                    (chord.lightcurve.emersion - chord.lightcurve.tref).sec + chord.lightcurve.exptime*40)
        pl.ylim(0 - chord.lightcurve.flux.std(ddof=1), np.median(chord.lightcurve.flux)+chord.lightcurve.flux.std(ddof=1)*3)
        pl.legend(loc=1, ncol=3)
        pl.tight_layout()
        pl.savefig(filename, dpi=400)
        pl.clf() 
        print("\tFigure '{}' was generated!".format(
            filename
        ))
    
def _projected_chords_plot(event):
        import matplotlib.pyplot as pl
        from sora.extra import draw_ellipse

        print("Projecting chords onto sky plane...")
    
        pl.figure(figsize=(8,4))
        event.chords.plot_chords(segment = 'positive', lw = 2, alpha=1, zorder=2)
        event.chords.plot_chords(segment = 'negative', lw = 2, alpha=0.5, zorder=2)
        event.chords.plot_chords(segment = 'error', lw = 4, alpha=1, zorder=1, color='r')

        # Mudar essa parte, aplicando o offset do centro do corpo...
        try:
            pl.xlim(-event.body.diameter.value, event.body.diameter.value)
            pl.ylim(-event.body.diameter.value, event.body.diameter.value)
        except:
            pl.xlim(-1500, 1500)
            pl.ylim(-1500, 1500)
    
        try:
            draw_ellipse(equatorial_radius=event.fitted_params['equatorial_radius'][0],
                         oblateness=event.fitted_params['oblateness'][0],
                         center_f =event.fitted_params['center_f'][0],
                         center_g=event.fitted_params['center_g'][0],
                         position_angle=event.fitted_params['position_angle'][0],
                         zorder=0, lw=1)
        except:
            pass
            
        pl.legend(bbox_to_anchor=(1.3,1.015))
        pl.tight_layout()    
        pl.savefig('Figures/projected_chords_plot.pdf', dpi=400)
        pl.clf()
        print("\tFigure 'Figures/projected_chords_plot.pdf' was generated!")