from byc import constants
from basico import *



output_selection = [
    'Time', 'Prot.ParticleNumber', 'YP.ParticleNumber',
    'Ydj1.ParticleNumber', 'ProtF.ParticleNumber', 'ProtM.ParticleNumber',
    'Oli.ParticleNumber', 'Agg.ParticleNumber', 'YOAgg.ParticleNumber',
    'YM.ParticleNumber', 'YO.ParticleNumber', 'Hsp104.ParticleNumber',
    'Cln3.ParticleNumber', 'YC.ParticleNumber', 'Cln3F.ParticleNumber',
    'Whi5.ParticleNumber', 'Whi5i.ParticleNumber',
    'Compartments[cell].Volume', 'Values[Budding]', 'Values[Generation]',
    'Values[TimeAtBudding]', 'Values[TimeNotinG1]',
    'Values[UnsyncedGeneration]', 'Values[Timer]'
]



def main():
    # Just get the model directly via its download URL on the biomodels website
    print('Loading model')
    model = load_model('https://www.ebi.ac.uk/biomodels/model/download/MODEL1901210001.5?filename=ProteinAggregationCellCycle.cps')
    output_dfs = []
    output_savepaths = []
    n_cells = 75
    for i in range(n_cells):
        print(f'Simulating cell {i+1} of {n_cells}')
        output = run_time_course_with_output(output_selection=output_selection,
                                    use_initial_values=True,
                                    method='directMethod',
                                    duration=1000,
                                    automatic=False,
                                    step_number=100)
        output.loc[:, 'cell_index'] = i
        savepath = os.path.join(constants.byc_data_dir, f'meta/ProteinAggregationCellCycle_cell{str(i).zfill(3)}_output.csv')
        output.to_csv(savepath)
        print(f'Saved model output at\n{savepath}')
        output_savepaths.append(savepath)
        output_dfs.append(output)

if __name__=='__main__':
    main()