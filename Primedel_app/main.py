from flask import Flask, flash, request, url_for,render_template,send_from_directory
from werkzeug.utils import secure_filename
import os,sys
import numpy as np
import primedel.design as pde



model=np.load('/'.join([os.path.dirname(pde.__file__),'indel_ratio.npz']))


app = Flask(__name__)
#CORS(app,supports_credentials=True)

ALLOWED_EXTENSIONS = {'txt', 'xls', 'scored'}
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[-1].lower() in ALLOWED_EXTENSIONS
app.secret_key = "secret key"
tmp_dir = '/tmp/'
#os.makedirs(tmp_dir, exist_ok = True)

@app.route("/", methods=['POST'])
def webtool():
    args = request.form.to_dict()
    args['pos_range'] = [int(s) for s in args['pos_range'].split()]
    args['size_range'] = [int(s) for s in args['size_range'].split()]
    args['precise'] = bool(args['precise'])
    args['homology_len'] = int(args['homology_len' ])
    args['threshold'] = int(args['threshold' ])
    seq = args['sequence'] 
    scored_file = request.files['scored_file']
    if scored_file and allowed_file(scored_file.filename):
        filename = secure_filename(scored_file.filename)
        if '.scored' in filename:
            guides = pde.read_flashfry(scored_file, args['threshold'])
        elif '.txt' in filename:
            guides = pde.read_gpp_designer(scored_file)
        elif '.xls' in filename:
            # Excel file is a bit weird-- not able to be read by pandas directly.
            fpath = os.path.join(tmp_dir, secure_filename(scored_file.filename))
            scored_file.save(fpath)
            guides = pde.read_crispor(fpath, args['threshold'])
    else:
        guides = pde.gen_guides(seq)
        #scored_file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
    fds = [(s[0], s[1],  # get guide and nick site
        13 + (pde.softmax(np.dot(pde.onehotencoder(s[0]),model['weights'])+model['bias'])[0]<0.75)//1 + # if ratio nick at 16 is >0.25
        ((s[0].upper().count('G')+s[0].upper().count('C'))/20<0.4)//1,s[-1])  # if the GC content is lower than 0.4 get PBS length
       for s in guides if s[2]=='FWD' and s[1]>30 ] 
    rvs = [(s[0],s[1],  # get guide and nick site
        13 + (pde.softmax(np.dot(pde.onehotencoder(s[0]),model['weights'])+model['bias'])[0]<0.75)//1 + 
        ((s[0].upper().count('G')+s[0].upper().count('C'))/20<0.4)//1,s[-1])
       for s in guides if s[2]=='REV' ]

    if args['pos_range']:
        if args['precise']:
            pairs = pde.peg_design_by_start_end(fds,rvs,seq,args['pos_range'],args['homology_len'],p=args['precise'])
        else:
            pairs = pde.peg_design_by_start_end(fds,rvs,seq,args['pos_range'],args['homology_len'])
    else:
        pairs = pde.peg_design_by_size(fds,rvs,seq,args['size_range'],args['homology_len'])
    # Writing stats files and fastq files 
    fname = os.path.join(tmp_dir, 'pegpairs_design.txt') 
    np.savetxt(fname,pairs,delimiter='\t',fmt='%s',
          header='RNA_1\thomology+PBS_1\tnick_1\tgRNA_2\thomology+PBS_2\tnick_2\tComments\tDeletion size\tExpected deletion result',comments='')
    return send_from_directory(tmp_dir, 'pegpairs_design.txt', as_attachment=True)

@app.route('/')
def static_page():
    return render_template('index.html')

if __name__ == "__main__":
    app.run(debug = True)