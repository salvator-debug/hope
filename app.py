from flask import Flask, render_template, request
import os
from flask import Flask, request, send_file
import pandas as pd
from io import BytesIO
from datetime import datetime
from flask import send_file, jsonify
import io
import json
import pdfkit
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle
from reportlab.lib import colors





app = Flask(__name__)

def calculate_reagents(summA, summT, summG, summC, summ, scale):
    if scale == '10мм':
        resA = round(summA * 100 / 19, 2)
        resT = round(summT * 100 / 26, 2)
        resG = round(summG * 100 / 20, 2)
        resC = round(summC * 100 / 22, 2)
        tetrazole = round(summ * 180 / 266, 2)
        tetrazole_percent = round(100 * tetrazole / 180, 2)
        anhydride = round(summ * 180 / 153, 2)
        anhydride_percent = round(100 * anhydride / 180, 2)
        nmi = round(summ * 180 / 165, 2)
        nmi_percent = round(100 * nmi / 180, 2)
        tca = round(summ * 450 / 106, 2)
        tca_percent = round(100 * tca / 450, 2)
        iodine = round(summ * 200 / 115, 2)
        iodine_percent = round(100 * iodine / 200, 2)
        acn = round(summ * 4000 / 130, 2)
        acn_percent = round(100 * acn / 4000, 2)
    elif scale == '1мм':
        resA = round(summA * 100 / 71, 2)
        resT = round(summT * 100 / 81, 2)
        resG = round(summG * 100 / 71, 2)
        resC = round(summC * 100 / 76, 2)
        tetrazole = round(summ * 180 / 417, 2)
        tetrazole_percent = round(100 * tetrazole / 180, 2)
        anhydride = round(summ * 180 / 621, 2)
        anhydride_percent = round(100 * anhydride / 180, 2)
        nmi = round(summ * 180 / 692, 2)
        nmi_percent = round(100 * nmi / 180, 2)
        tca = round(summ * 450 / 385, 2)
        tca_percent = round(100 * tca / 450, 2)
        iodine = round(summ * 200 / 488, 2)
        iodine_percent = round(100 * iodine / 200, 2)
        acn = round(summ * 4000 / 500, 2)
        acn_percent = round(100 * acn / 4000, 2)
    elif scale == '0.2мм':
        resA = round(summA * 100 / 107, 2)
        resT = round(summT * 100 / 124, 2)
        resG = round(summG * 100 / 107, 2)
        resC = round(summC * 100 / 104, 2)
        tetrazole = round(summ * 180 / 446, 2)
        tetrazole_percent = round(100 * tetrazole / 180, 2)
        anhydride = round(summ * 180 / 621, 2)
        anhydride_percent = round(100 * anhydride / 180, 2)
        nmi = round(summ * 180 / 692, 2)
        nmi_percent = round(100 * nmi / 180, 2)
        tca = round(summ * 450 / 385, 2)
        tca_percent = round(100 * tca / 450, 2)
        iodine = round(summ * 200 / 488, 2)
        iodine_percent = round(100 * iodine / 200, 2)
        acn = round(summ * 4000 / 500, 2)
        acn_percent = round(100 * acn / 4000, 2)
    elif scale == 'праймер':
        resA = round(summA * 100 / 214, 2)
        resT = round(summT * 100 / 248, 2)
        resG = round(summG * 100 / 214, 2)
        resC = round(summC * 100 / 208, 2)
        tetrazole = round(summ * 180 / 750, 2)
        tetrazole_percent = round(100 * tetrazole / 180, 2)
        anhydride = round(summ * 180 / 621, 2)
        anhydride_percent = round(100 * anhydride / 180, 2)
        nmi = round(summ * 180 / 692, 2)
        nmi_percent = round(100 * nmi / 180, 2)
        tca = round(summ * 450 / 385, 2)
        tca_percent = round(100 * tca / 450, 2)
        iodine = round(summ * 200 / 488, 2)
        iodine_percent = round(100 * iodine / 200, 2)
        acn = round(summ * 4000 / 570, 2)
        acn_percent = round(100 * acn / 4000, 2)
    else:
        raise ValueError(f"Невідомий масштаб: {scale}")



    return {
        'resA': resA,
        'resT': resT,
        'resG': resG,
        'resC': resC,
        'tetrazole': tetrazole,
        'tetrazole_percent': tetrazole_percent,
        'anhydride': anhydride,
        'anhydride_percent': anhydride_percent,
        'nmi': nmi,
        'nmi_percent': nmi_percent,
        'tca': tca,
        'tca_percent': tca_percent,
        'iodine': iodine,
        'iodine_percent': iodine_percent,
        'acn': acn,
        'acn_percent': acn_percent
    }

@app.route('/', methods=['GET', 'POST'])
def index():
    results = []
    total = {
        'A': 0, 'T': 0, 'G': 0, 'C': 0,
        'tetrazole': 0, 'tetrazole_percent': 0,
        'anhydride': 0, 'anhydride_percent': 0,
        'nmi': 0, 'nmi_percent': 0,
        'tca': 0, 'tca_percent': 0,
        'iodine': 0, 'iodine_percent': 0,
        'acn': 0, 'acn_percent': 0,
        'summmolarweight': 0
    }

    if request.method == 'POST':
        scale = request.form.get('scale')
        for i in range(1, 5):
            if request.form.get(f'use_seq{i}'):
                seq_raw = request.form.get(f'sequence{i}', '').upper().replace(" ", "").replace("\n", "")
                summA = seq_raw.count('A')
                summT = seq_raw.count('T')
                summG = seq_raw.count('G')
                summC = seq_raw.count('C')
                summ = summA + summT + summG + summC
                summmolarweight = summA * 313.2 + summT * 304.2 + summG * 329.2 + summC * 289.2

                reagents = calculate_reagents(summA, summT, summG, summC, summ, scale)
                if not reagents:
                    continue

                results.append({
                    'seq_num': i,
                    'sequence': seq_raw,
                    'A': reagents['resA'],
                    'T': reagents['resT'],
                    'G': reagents['resG'],
                    'C': reagents['resC'],
                    'tetrazole': reagents['tetrazole'],
                    'tetrazole_percent': reagents['tetrazole_percent'],
                    'anhydride': reagents['anhydride'],
                    'anhydride_percent': reagents['anhydride_percent'],
                    'nmi': reagents['nmi'],
                    'nmi_percent': reagents['nmi_percent'],
                    'tca': reagents['tca'],
                    'tca_percent': reagents['tca_percent'],
                    'iodine': reagents['iodine'],
                    'iodine_percent': reagents['iodine_percent'],
                    'acn': reagents['acn'],
                    'acn_percent': reagents['acn_percent'],
                    'summmolarweight': round(summmolarweight, 2)
                })

                # Сумуємо до total
                total['A'] += reagents['resA']
                total['T'] += reagents['resT']
                total['G'] += reagents['resG']
                total['C'] += reagents['resC']
                total['tetrazole'] += reagents['tetrazole']
                total['tetrazole_percent'] += reagents['tetrazole_percent']
                total['anhydride'] += reagents['anhydride']
                total['anhydride_percent'] += reagents['anhydride_percent']
                total['nmi'] += reagents['nmi']
                total['nmi_percent'] += reagents['nmi_percent']
                total['tca'] += reagents['tca']
                total['tca_percent'] += reagents['tca_percent']
                total['iodine'] += reagents['iodine']
                total['iodine_percent'] += reagents['iodine_percent']
                total['acn'] += reagents['acn']
                total['acn_percent'] += reagents['acn_percent']
                total['summmolarweight'] += round(summmolarweight, 2)

        # Округлюємо total
        for key in total:
            total[key] = round(total[key], 2)

    return render_template('index.html', results=results, total=total)











# ___________________________________________________________________________________________________________

@app.route('/edit', methods=['GET', 'POST'])
def edit():
    count = 1 #кількість рядків у таблиці
    columns = [
        'Дата', 'Аргон', 'Ацетонітрил', 'Температура',
        'G', 'A', 'C', 'T', '6', '7', '8', 'Тетразол',
        '10', 'Кислотний ангідрид', 'NMI', '13', 'TCA', 'Суміш'
    ]

    initial_value_raw = request.form.get('initial_value')
    initial_value = None
    try:
        if initial_value_raw is not None:
            initial_value = float(initial_value_raw)
    except ValueError:
        initial_value = None

    # Індивідуальні значення для першого рядка (по колонках)
    initial_row_values = {}
    for col in columns:
        val = request.form.get(f'initial_{col}')
        if val is not None and val != '':
            initial_row_values[col] = val

    # Віднімання для кожної колонки (по дефолту 10)
    subtract_values = {}
    for col in columns:
        val = request.form.get(f'subtract_{col}', '10')
        try:
            subtract_values[col] = float(val)
        except ValueError:
            subtract_values[col] = 10

    results = []

    for i in range(count):
        row = {}
        for col in columns:
            key = f"{col}_{i}"
            if i == 0:
                if col != 'Дата':
                    if col in initial_row_values:
                        try:
                            v = float(initial_row_values[col])
                            row[col] = str(v - subtract_values.get(col, 10))
                        except ValueError:
                            row[col] = initial_row_values[col]
                    elif initial_value is not None:
                        row[col] = str(initial_value - subtract_values.get(col, 10))
                    else:
                        row[col] = request.form.get(key, '100')
                else:
                    row[col] = request.form.get(key, datetime.today().strftime('%Y-%m-%d'))
            else:
                row[col] = request.form.get(key, '')
        results.append(row)

    return render_template('edit.html',
                           results=results,
                           columns=columns,
                           initial_value=initial_value_raw or '100',
                           initial_row_values=initial_row_values,
                           subtract_values=subtract_values)

@app.route('/export/excel', methods=['POST'])
def export_excel():
    count = int(request.form.get('count', 0))
    data = []
    for i in range(count):
        row = {
            'Дата': request.form.get(f'date_{i}', ''),
            'Аргон': '-',
            'Ацетонітрил (мл)': request.form.get(f'acn_{i}', ''),
            'Температура': '-',
            'G %': request.form.get(f'G_{i}', ''),
            'A %': request.form.get(f'A_{i}', ''),
            'C %': request.form.get(f'C_{i}', ''),
            'T %': request.form.get(f'T_{i}', ''),
            '6': '-',
            '7': '-',
            '8': '-',
            'Тетразол (мл)': request.form.get(f'tetrazole_{i}', ''),
            'Кислотний ангідрид (мл)': request.form.get(f'anhydride_{i}', ''),
            'NMI (мл)': request.form.get(f'nmi_{i}', ''),
            '13': '-',
            'TCA (мл)': request.form.get(f'tca_{i}', ''),
            'Суміш': '-'
        }
        data.append(row)
    df = pd.DataFrame(data)
    output = BytesIO()
    df.to_excel(output, index=False)
    output.seek(0)

    filename = f"results_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.xlsx"
    return send_file(output, as_attachment=True, download_name=filename, mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
@app.route('/export/pdf', methods=['POST'])
def export_pdf():
    count = int(request.form.get('count', 0))
    data = []
    headers = ['Дата', 'Аргон', 'Ацетонітрил', 'Температура', 'G', 'A', 'C', 'T',
               '6', '7', '8', 'Тетразол', 'Кислотний ангідрид', 'NMI', '13', 'TCA', 'Суміш']
    data.append(headers)
    for i in range(count):
        row = [
            request.form.get(f'date_{i}', ''),
            '-',
            request.form.get(f'acn_{i}', ''),
            '-',
            request.form.get(f'G_{i}', ''),
            request.form.get(f'A_{i}', ''),
            request.form.get(f'C_{i}', ''),
            request.form.get(f'T_{i}', ''),
            '-',
            '-',
            '-',
            request.form.get(f'tetrazole_{i}', ''),
            request.form.get(f'anhydride_{i}', ''),
            request.form.get(f'nmi_{i}', ''),
            '-',
            request.form.get(f'tca_{i}', ''),
            '-'
        ]
        data.append(row)
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    table = Table(data, repeatRows=1)
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#f0f0f0')),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.black),
        ('FONTSIZE', (0, 0), (-1, -1), 8),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER')
    ])
    table.setStyle(style)
    doc.build([table])
    buffer.seek(0)
    filename = f"results_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.pdf"
    return send_file(buffer, as_attachment=True, download_name=filename, mimetype='application/pdf')
# ___________________________________________________________________________________________________________
@app.route('/mode2', methods=['GET', 'POST'])
def mode2():
    total = {
        'A': 0, 'T': 0, 'G': 0, 'C': 0,
        'tetrazole': 0,'tetrazole_percent':0,
        'anhydride_percent':0,'nmi_percent':0,
        'tca_percent':0,'iodine_percent':0,
        'acn_percent':0, 'anhydride': 0,
        'nmi': 0, 'tca': 0, 'iodine': 0,
        'acn': 0, 'summmolarweight': 0
    }
    if request.method == 'POST':
        sequences = request.form.getlist('sequence')
        scales = request.form.getlist('scale')
        for seq, scale in zip(sequences, scales):
            seq = seq.upper().replace(" ", "").replace("\n", "")
            summA = seq.count('A')
            summT = seq.count('T')
            summG = seq.count('G')
            summC = seq.count('C')
            summ = summA + summT + summG + summC
            summmolarweight = summA * 313.2 + summT * 304.2 + summG * 329.2 + summC * 289.2
            reagents = calculate_reagents(summA, summT, summG, summC, summ, scale)
            if reagents:
                total['A'] = round(total['A']+reagents['resA'],2)
                total['T'] = round(total['T']+reagents['resT'],2)
                total['G'] = round(total['G']+reagents['resG'],2)
                total['C'] = round(total['C']+reagents['resC'],2)
                total['tetrazole'] = round(total['tetrazole'] + reagents['tetrazole'],2)
                total['tetrazole_percent'] = round(total['tetrazole_percent'] + reagents['tetrazole_percent'],2)
                total['anhydride'] = round(total['anhydride'] + reagents['anhydride'], 2)
                total['anhydride_percent'] = round(total['anhydride_percent'] + reagents['anhydride_percent'], 2)
                total["nmi"] = round(total['nmi']+reagents['nmi'],2)
                total['nmi_percent'] = round(total['nmi_percent'] + reagents['nmi_percent'],2)
                total['tca'] = round(total['tca']+reagents['tca'],2)
                total['tca_percent']  = round(total['tca_percent'] + reagents['tca_percent'],2)
                total['iodine'] = round(total['iodine']+reagents['iodine'],2)
                total['iodine_percent'] = round(total['iodine_percent']+reagents['iodine_percent'],2)
                total['acn'] = round(total['acn']+reagents['acn'],2)
                total['acn_percent'] = round(total['acn_percent']+reagents['acn_percent'],2)
                total['summmolarweight'] = round(total['summmolarweight']+summmolarweight,2)
    return render_template('mode2_content.html', total=total)
# ___________________________________________________________________________________________________________

from flask import Flask, render_template, request



# Енергії для пар нуклеотидів (ккал/моль)
# Значення приблизно взяті з nearest-neighbor моделі (SantaLucia, 1998)
NN_ENERGIES = {
    ("A", "T"): -1.0,
    ("T", "A"): -1.0,
    ("G", "C"): -1.5,
    ("C", "G"): -1.5
}

def reverse_complement(seq):
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement.get(base, "") for base in reversed(seq))

def calc_delta_g(seq1, seq2):
    # Прибираємо пробіли, переводи рядків і переводимо в верхній регістр
    seq1 = seq1.replace(" ", "").replace("\n", "").strip().upper()
    seq2 = seq2.replace(" ", "").replace("\n", "").strip().upper()
    valid_bases = set("ATGC")

    # Перевірка валідності
    if not seq1 or not seq2:
        return 0.0
    if len(seq1) < 4 or len(seq2) < 4:
        return 0.0
    if not (set(seq1) <= valid_bases and set(seq2) <= valid_bases):
        return 0.0

    # Реверсно-комплементарна послідовність
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rc_seq2 = ''.join(complement[b] for b in reversed(seq2))

    # Пошук максимальної кількості пар
    max_pairs = 0
    for offset in range(-len(rc_seq2), len(seq1)):
        pairs = 0
        for i in range(len(seq1)):
            j = i - offset
            if 0 <= j < len(rc_seq2):
                if seq1[i] == rc_seq2[j]:
                    pairs += 1
        max_pairs = max(max_pairs, pairs)

    # Формула ΔG
    delta_g = -1.5 * max_pairs
    return round(delta_g, 2)

def calculate_highest_delta_g(systems):
    results = []
    for system in systems:
        dna_list = system['dna']
        min_dg = None
        pairs = []
        n = len(dna_list)
        for i in range(n):
            for j in range(i + 1, n):
                dg = calc_delta_g(dna_list[i], dna_list[j])
                if dg is not None:
                    pairs.append((dna_list[i], dna_list[j], dg))
                    if min_dg is None or dg < min_dg:
                        min_dg = dg
        results.append({
            'system_name': system['name'],
            'max_delta_g': min_dg,
            'pairs': pairs
        })
    return results

@app.route('/mode3', methods=['GET', 'POST'])
def mode3():
    systems = []
    results = []

    if request.method == 'POST':
        form = request.form
        system_indices = sorted(set(k.split('_')[1] for k in form if k.startswith('virus_')))
        for idx in system_indices:
            name = form.get(f'virus_{idx}', f'Система {idx}')
            dna_seqs = form.getlist(f'dna_{idx}[]')
            dna_seqs = [seq.strip().upper() for seq in dna_seqs if seq.strip()]
            systems.append({'name': name, 'dna': dna_seqs})

        results = calculate_highest_delta_g(systems)

    return render_template('mode3_content.html', systems=systems, results=results)

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port)