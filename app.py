from flask import Flask, render_template, request

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
from flask import send_file, jsonify
import pandas as pd
import io
import json

@app.route('/edit', methods=['POST'])
def edit():
    results_json = request.form.get('results_json')
    if not results_json:
        return "No data provided", 400
    results = json.loads(results_json)
    return render_template('edit.html', results=results)

@app.route('/export', methods=['POST'])
def export():
    count = int(request.form.get('count', 0))
    data = []
    for i in range(count):
        row = {
            "Колонка": request.form.get(f'seq_num_{i}'),
            "Послідовність": request.form.get(f'sequence_{i}'),
            "A %": request.form.get(f'A_{i}'),
            "T %": request.form.get(f'T_{i}'),
            "G %": request.form.get(f'G_{i}'),
            "C %": request.form.get(f'C_{i}'),
            "Тетразол (мл)": request.form.get(f'tetrazole_{i}'),
            "Кислотний ангідрид (мл)": request.form.get(f'anhydride_{i}'),
            "NMI (мл)": request.form.get(f'nmi_{i}'),
            "TCA (мл)": request.form.get(f'tca_{i}'),
            "Йод (мл)": request.form.get(f'iodine_{i}'),
            "Ацетонітрил (мл)": request.form.get(f'acn_{i}'),
            "Молярна маса (г/моль)": request.form.get(f'summmolarweight_{i}')
        }
        data.append(row)

    df = pd.DataFrame(data)
    output = io.BytesIO()
    df.to_excel(output, index=False, engine='openpyxl')
    output.seek(0)

    return send_file(output, as_attachment=True, download_name="results.xlsx", mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')

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

if __name__ == '__main__':
    app.run(debug=True)