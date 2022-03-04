"""
Automate ALE analyses to integrate over a distribution of
host species trees. Can use either the dated or undated
version of ALE, and can evaluate the DTL model and all subsets 
of the DTL model when model choice is requested.   
author: J. Satler
date: 2 Feb 2019
version: 3
usage:
    python autoALE.py figTrees WaspTrees traits dated|undated model choice? [y|n]
"""
import os
import sys
import shutil
import dendropy
import datetime
import subprocess
import pandas as pd

models = [["DTL"],
          ["DL", "tau=0"],
          ["TL", "delta=0"],
          ["DT", "lambda=0"],
          ["L", "delta=0", "tau=0"],
          ["D", "tau=0", "lambda=0"],
          ["T", "delta=0", "lambda=0"],
          ["CoSp", "delta=0", "tau=0", "lambda=0"]]

def read_beast_trees(tr):
    """read a set of newick rooted species trees"""
    schema = 'newick'

    #is nexus or newick tree file?
    with open(tr, 'r') as t:
        if '#NEXUS' in t.readline():
            schema = 'nexus'
    return dendropy.TreeList.get(file=open(tr, 'r'),
                                 schema=schema,
                                 tree_offset=0,
                                 rooting="default-rooted")

def burnin_and_thin_dist(trees):
    """burnin 10 percent of trees and thin to 1000 total trees"""
    #burnin distribution
    bi_trees = trees[int(len(trees) * .10):]

    #return all post-burnin trees if distribution fewer than 1000 trees
    if len(bi_trees) <= 1000:
        return bi_trees

    new_trees = []
    tree = 0
    for i in range(len(bi_trees)):
        if i == int(round(tree)):
            new_trees.append(bi_trees[i])
            #thin distribution to 1000 trees
            tree += len(bi_trees) / 1000.0
    return new_trees

def write_beast_tree(tree):
    """write a species tree to file for analysis"""
    with open("SpeciesTree.tre", 'w') as out:
        out.write(tree.as_string('newick'))

def get_associations(traits):
    """get associations of gene tree tips to species tree tips"""
    with open(traits, 'r') as tr:
        return {t.split()[0].strip():t.split()[1].strip() for t in tr
                                                          if not t.startswith("#")}

def assign_tips(file, traits):
    """assign gene tree tips to species tree tips"""
    new_trees = []
    with open(file, 'r') as t:
        for line in t:
            #replace tip names with associations
            for k, v in traits.items():
                line = line.replace(k, v + '_' + k)
            new_trees.append(line)
    #write new trees to file
    with open(file.split("/")[-1].split(".")[0] + "_READY.trees", 'w') as out:
        for i in new_trees:
            out.write(i)
    return file.split("/")[-1].split(".")[0] + "_READY.trees"

def run_ale_observe(tr):
    """run ALEobserve"""
    # cmd = "docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve " + tr + " burnin=1000"
    cmd = ["docker", "run", "-v", os.getcwd() + ":" + os.getcwd(), "-w", os.getcwd(), "boussau/alesuite", "ALEobserve", tr, "burnin=3000"]    # subprocess.call(cmd, shell=True)
    subprocess.call(cmd)

def run_ale(genetrees, ale, model):
    """run ALEml"""
    if ale == "undated":
        cmd = ["docker", "run", "-v", os.getcwd() + ":" + os.getcwd(), "-w", os.getcwd(), "boussau/alesuite", "ALEml_undated" , "SpeciesTree.tre", genetrees + '.ale', "sample=1000"]
        # cmd = "docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated SpeciesTree.tre " + genetrees + ".ale sample=1000"
        if len(model) > 1:
            cmd.extend(model[1:])
            # cmd = cmd + model[1:]
    else:
        cmd = ["docker", "run", "-v", os.getcwd() + ":" + os.getcwd(), "-w", os.getcwd(), "boussau/alesuite", "ALEml" , "SpeciesTree.tre", genetrees + '.ale', "sample=1000"]
        # cmd = "docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml SpeciesTree.tre " + genetrees + ".ale sample=1000"
        if len(model) > 1:
            cmd.extend(model[1:])
            # cmd = cmd + model[1:]

    #run ale
    # p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    #retain log likelihood in case program exits with error
    out, err = p.communicate()
    for line in out.split("\n"):
        if line.startswith("LL="):
            return line[3:]
    
    # TMP output err
    # if err:
    #    print err

def get_results(folder):
    """get parameter estimates from ALE"""
    output = [f for f in os.listdir('.') if f.endswith('ml_rec')]

    #return if ALE threw an error
    if not output:
        return

    with open(output[0], 'r') as res:
        stats = {}
        for line in res:
            if line.startswith(">logl"):
                stats['logl'] = line.strip().split( )[1]
            elif line.startswith("ML"):
                line = line.strip().split( )
                stats['MLrates'] = line[1:]
            elif line.startswith("Total"):
                line = line.strip().split( )
                stats['Total'] = line[1:]
        return stats

def set_file_dir(folder, model, num):
    """create and place run files in directory"""
    cur = os.path.realpath('.').split('/')[:-1]
    out = '/'.join(cur) + '/results/output_' + folder + '/' + model + '/tree' + str(num)

    #create and populate output folders for each species tree
    os.makedirs(out)
    ALEfiles = [shutil.move(f, out) for f in os.listdir('.')
                                          if '.ale.' in f]
    GTfiles = [shutil.copy(f, out) for f in os.listdir('.')
                                         if '_READY.trees' in f]
    shutil.copy('SpeciesTree.tre', out)

def summarize_mc(resultsD, folder):
    """place results in output file"""
    modsel = {}
    for k, v in resultsD.items():
        if v:
            with open("results_summary_" + k + ".txt", 'w') as out:
                header = "SpTree\tlogl\tdelta\ttau\tlambda\tduplications\ttransfers\tlosses\tspeciations"
                out.write(header + '\n')
                modsel[k] = [header.split('\t')]

                for i in v:
                    if len(i) == 2:
                        res = 'tree' + i["tree"] + '\t' + i['logl'] + "\t0\t0\t0\t0\t0\t0\t0"
                    else:
                        res = 'tree' + i["tree"] + '\t' + i['logl'] + '\t' + '\t'.join(i['MLrates']) + '\t' + '\t'.join(i['Total'])
                    out.write(res + '\n')
                    modsel[k].append(res.split('\t'))
            #place summary file in respective folder
            if os.path.isdir('../results/output_' + folder + '/' + k):
                shutil.move("results_summary_" + k + ".txt", '../results/output_' + folder + '/' + k)
            else:
                os.makedirs('../results/output_' + folder + '/' + k)
                shutil.move("results_summary_" + k + ".txt", '../results/output_' + folder + '/' + k)
    return modsel

def summarize_modsel(models, folder):
    """summarize parameters for model selection"""
    modsel_table = []
    for k, v in sorted(models.items()):
        data = pd.DataFrame(v[1:], columns=v[0], dtype=float)

        #get header
        if not modsel_table:
            modsel_table.append(("model", data.columns.values[1:].tolist()))

        vals = [str(content.mean()) for label, content in data.iteritems() if label != "SpTree"]
        modsel_table.append((k, vals))

    #write to file
    with open("_model_choice_results.txt", "w") as out:
        for r in modsel_table:
            out.write(r[0] + "\t" + "\t".join(r[1]) + "\n")
    shutil.move("_model_choice_results.txt", '../results/output_' + folder)

def main():
    if len(sys.argv) < 5:
        print "python autoALE.py figTrees figwaspTrees traits Dated/Undated model choice [y|n] burnins"
        sys.exit()

    #read BEAST species trees
    t = read_beast_trees(sys.argv[1])

    #burnin and thin distribution
    tr = burnin_and_thin_dist(t)

    #assign tips and get CCPs for fig wasp trees
    traits = get_associations(sys.argv[3])
    genetrees = assign_tips(sys.argv[2], traits)
    run_ale_observe(genetrees)

    #specify dated (default) or undated analysis
    ale = 'dated'
    if sys.argv[4].lower().startswith("u"):
        print 'undated analysis'
        ale = 'undated'
    else:
        print 'dated analysis'

    #total results
    resD = {m[0]:[] for m in models}

    for i in range(len(tr)):
        write_beast_tree(tr[i])
        print "{0}\tTree {1}".format(datetime.datetime.now(), i)

        for m in models:
            logL = run_ale(genetrees, ale, m)

            #write results
            results = get_results(ale)

            #if ALE threw an error
            if not results:
                #print "Rep {0} threw an ERROR and was SKIPPED".format(i)
                resD[m[0]].append({'tree':str(i), 'logl':logL})
                continue

            results['tree'] = str(i)
            resD[m[0]].append(results)

            set_file_dir(ale, m[0], i)

            #check if model choice is wanted
            #default is no unless specified
            if not sys.argv[5].lower().startswith("y"):
                break

        os.remove('SpeciesTree.tre')
        continue

    #write stats to file
    modsel = summarize_mc(resD, ale)
    #summarize model choice information
    summarize_modsel(modsel, ale)

    #clean directory
    rmfiles = [os.remove(f) for f in os.listdir('.') if '_READY.trees' in f]

if __name__ == '__main__':
    main()

