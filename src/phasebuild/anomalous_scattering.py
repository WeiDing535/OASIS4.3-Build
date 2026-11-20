"""Offline theoretical anomalous factors from explicitly supplied wavelength."""
import math
import re
import gemmi

SOURCE_URL = 'https://gemmi.readthedocs.io/en/stable/scattering.html#anomalous-scattering'


def calculate(elements, wavelength):
    if isinstance(wavelength, bool) or not isinstance(wavelength, (int, float)) or not math.isfinite(wavelength) or not .1 <= wavelength <= 10:
        raise ValueError('please provide the actual experimental SAD wavelength (0.1–10 Å); the wavelength cannot be left empty and will not be estimated automatically from the MTZ')
    if not isinstance(elements, str) or not elements.strip():
        raise ValueError('please provide heavy-atom elements, e.g. Se or S; separate multiple elements with commas')
    names = []
    for text in re.split(r'[,;\s]+', elements.replace('\uff0c', ',').strip()):
        element = gemmi.Element(text)
        if not 3 <= element.atomic_number <= 92 or not re.fullmatch('[A-Za-z]{1,2}', text):
            raise ValueError('unsupported heavy-atom element: ' + text + '; theoretical calculation supports atomic numbers 3–92')
        if element.atomic_number >= 79 and wavelength < .16:
            raise ValueError('this heavy element does not fit the current theoretical model in this short wavelength range; please check the experimental wavelength')
        if element.name not in names:
            names.append(element.name)
    energy = gemmi.hc / wavelength
    factors = {}
    for name in names:
        fp, fpp = gemmi.cromer_liberman(gemmi.Element(name).atomic_number, energy)
        if not math.isfinite(fp) or not math.isfinite(fpp) or fpp <= 0:
            raise ValueError('no valid theoretical f″ can be obtained for this element and wavelength: ' + name)
        factors[name] = {'fp': fp, 'fpp': fpp}
    return {'mode': 'automatic_theoretical', 'elements': names, 'wavelength_angstrom': wavelength,
            'energy_eV': energy, 'factors': factors,
            'fpp_argument': ','.join(name + '=' + repr(factors[name]['fpp']) for name in names),
            'source': 'Gemmi Cromer–Liberman (Kissel–Pratt correction)',
            'gemmi_version': gemmi.__version__, 'source_url': SOURCE_URL,
            'wavelength_source': 'explicit user input; never inferred from MTZ',
            'scope': 'Theoretical isolated-atom values for the P+ substructure calculation. '
                     'Phaser computes its SAD factors internally from the supplied wavelength.',
            'limitation': 'Theoretical values are not an experimental absorption-edge scan.'}
