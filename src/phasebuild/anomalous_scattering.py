"""Offline theoretical anomalous factors from explicitly supplied wavelength."""
import math
import re
import gemmi

SOURCE_URL = 'https://gemmi.readthedocs.io/en/stable/scattering.html#anomalous-scattering'


def calculate(elements, wavelength):
    if isinstance(wavelength, bool) or not isinstance(wavelength, (int, float)) or not math.isfinite(wavelength) or not .1 <= wavelength <= 10:
        raise ValueError('请填写实际实验 SAD 波长（0.1–10 Å）；波长不能留空，也不会从 MTZ 自动估算')
    if not isinstance(elements, str) or not elements.strip():
        raise ValueError('请填写重原子元素，例如 Se 或 S，多元素用逗号分隔')
    names = []
    for text in re.split(r'[,;\s]+', elements.replace('，', ',').strip()):
        element = gemmi.Element(text)
        if not 3 <= element.atomic_number <= 92 or not re.fullmatch('[A-Za-z]{1,2}', text):
            raise ValueError('不支持的重原子元素：' + text + '；自动理论计算支持原子序数 3–92')
        if element.atomic_number >= 79 and wavelength < .16:
            raise ValueError('该重元素在此短波长范围不适用当前理论模型，请检查实验波长')
        if element.name not in names:
            names.append(element.name)
    energy = gemmi.hc / wavelength
    factors = {}
    for name in names:
        fp, fpp = gemmi.cromer_liberman(gemmi.Element(name).atomic_number, energy)
        if not math.isfinite(fp) or not math.isfinite(fpp) or fpp <= 0:
            raise ValueError('该元素与波长无法得到有效的理论 f″：' + name)
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
