import sys
import os

# need 2 file: mer_counts.jf and gtouch (also k)
k = int(sys.argv[1])
GTOUCH_FILE = sys.argv[2]
INVALID_UNITIG_FILE=""
GENOME_FILE=""
TOUCHED_GENOME_ID_FILE=""
if(len(sys.argv)>3):
    TOUCHED_GENOME_ID_FILE = GTOUCH_FILE
    INVALID_UNITIG_FILE = sys.argv[3]
    GENOME_FILE = sys.argv[4]

# ref = sys.argv[3]
colors = ["black", "gray", "aqua", "blue", "blueviolet", "brown", "burlywood", "cadetblue", "chartreuse", "chocolate",
          "coral", "cornflowerblue", "cornsilk", "crimson", "darkblue", "darkorchid", "aquamarine", "cyan", "darkcyan",
          "darkgoldenrod", "darkgreen", "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange", "darkred",
          "darksalmon", "darkseagreen", "darkslateblue", "darkturquoise", "darkviolet", "deeppink", "deepskyblue",
          "dodgerblue", "firebrick", "forestgreen", "fuchsia", "gainsboro", "gold", "goldenrod", "green", "greenyellow",
          "hotpink", "indianred", "indigo", "khaki", "lavender", "lavenderblush", "lawngreen", "lemonchiffon",
          "lightblue", "lightcoral", "lightcyan", "lightgoldenrodyellow", "lightpink", "lightsalmon", "lightseagreen",
          "lightskyblue", "lightsteelblue", "lightyellow", "limegreen", "linen", "magenta", "maroon",
          "mediumaquamarine", "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen", "mediumslateblue",
          "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue", "mintcream", "mistyrose",
          "moccasin", "navajowhite", "navy", "oldlace", "olive", "olivedrab", "orange", "orangered", "orchid",
          "palegoldenrod", "palegreen", "paleturquoise", "palevioletred", "papayawhip", "peachpuff", "peru", "pink",
          "plum", "powderblue", "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown", "seagreen", "sienna",
          "silver", "skyblue", "slateblue", "springgreen", "steelblue", "tan", "teal", "thistle", "tomato", "turquoise",
          "violet", "wheat", "yellow", "yellowgreen", "darkgray", "darkslategray", "dimgray", "lightslategray",
          "lightgray", "slategray"]




def read_gtouch(GTOUCH_FILE, kmer_to_id, id_to_kmer, genome_to_gid, genomes, genomefirst, genomelast, firsts, lasts, ufirsts, ulasts, unitigs):
    curr_id = 0
    genomeid = 0
    with open(GTOUCH_FILE, 'r') as reader:
        line = reader.readline()
        old_inv_uni = ""
        flag = False

        while line != '':
            line = line.rstrip()
            linesplit = line.split()
            inv_uni = linesplit[0].rstrip()
            if (inv_uni != old_inv_uni):
                unitigs.append(inv_uni)

                flag = True
            if (flag):
                for i in range(len(inv_uni) - k + 1):
                    kmer_to_append = inv_uni[i:i + k]
                    if not kmer_to_append in kmer_to_id:
                        kmer_to_id[kmer_to_append] = curr_id
                        id_to_kmer[curr_id] = kmer_to_append
                        curr_id += 1
            old_inv_uni = inv_uni

            g = linesplit[1].rstrip()

            # populating dict genome_to_gid
            if g not in genome_to_gid:
                genome_to_gid[g] = genomeid
                genomeid += 1
                genomes.append(g)

            firstkmer = inv_uni[0:k]
            lastkmer = inv_uni[len(inv_uni) - k:]

            for i in range(len(g) - k + 1):
                kmer_to_append = g[i:i + k]
                if not kmer_to_append in kmer_to_id:
                    kmer_to_id[kmer_to_append] = curr_id
                    id_to_kmer[curr_id] = kmer_to_append
                    curr_id += 1

                if (kmer_to_id[kmer_to_append] not in vertex_to_genomeid):
                    vertex_to_genomeid[kmer_to_id[kmer_to_append]] = set()
                vertex_to_genomeid[kmer_to_id[kmer_to_append]].add(genomeid)

            genomefirst[kmer_to_id[g[0:k]]] = genomeid
            genomelast[kmer_to_id[g[len(g) - k:len(g)]]] = genomeid

            firsts.append(kmer_to_id[firstkmer])
            lasts.append(kmer_to_id[lastkmer])

            line = reader.readline()

    for u in unitigs:
        ufirsts.add(kmer_to_id[u[0:k]])
        ulasts.add(kmer_to_id[u[len(u) - k:]])
    return curr_id



def read_gtouch_idversion(INVALID_UNITIG_FILE, GENOME_FILE, TOUCHED_GENOME_ID_FILE, kmer_to_id, id_to_kmer, genome_to_gid, genomes, genomefirst, genomelast, firsts, lasts, ufirsts, ulasts, unitigs):
    #feed invalid unitigs
    curr_id=0
    genomeid_reduced=0
    genomeid=-1

    reader_id=open(TOUCHED_GENOME_ID_FILE, 'r')

    with open(GENOME_FILE, 'r') as reader:
        line = reader.readline()
        id_of_touched_g = int(reader_id.readline().rstrip())
        while line != '':
            g = line.rstrip()
            if(g[0]!='>'):
                genomeid+=1
                if(id_of_touched_g == genomeid):
                    if g not in genome_to_gid:
                        genome_to_gid[g] = genomeid_reduced
                        genomes.append(g)
                    for i in range(len(g) - k + 1):
                        kmer_to_append = g[i:i + k]
                        if not kmer_to_append in kmer_to_id:
                            kmer_to_id[kmer_to_append] = curr_id
                            id_to_kmer[curr_id] = kmer_to_append
                            curr_id += 1

                        if (kmer_to_id[kmer_to_append] not in vertex_to_genomeid):
                            vertex_to_genomeid[kmer_to_id[kmer_to_append]] = set()
                        vertex_to_genomeid[kmer_to_id[kmer_to_append]].add(genomeid_reduced)

                    genomefirst[kmer_to_id[g[0:k]]] = genomeid_reduced
                    genomelast[kmer_to_id[g[len(g) - k:len(g)]]] = genomeid_reduced


                    genomeid_reduced+=1

                    id_of_touched_g = (reader_id.readline().rstrip())
                    if(id_of_touched_g==''):
                        break
            line = reader.readline()
    reader_id.close()
    with open(INVALID_UNITIG_FILE, 'r') as reader:
        line = reader.readline()
        while line != '':
            u = line.rstrip()
            if (u[0] != '>'):
                unitigs.append(u)
                for i in range(len(u) - k + 1):
                    kmer_to_append = u[i:i + k]
                    if not kmer_to_append in kmer_to_id:
                        kmer_to_id[kmer_to_append] = curr_id
                        id_to_kmer[curr_id] = kmer_to_append
                        curr_id += 1
                firstkmer = u[0:k]
                lastkmer = u[len(u) - k:]
                firsts.append(kmer_to_id[firstkmer])
                lasts.append(kmer_to_id[lastkmer])
                ufirsts.add(kmer_to_id[u[0:k]])
                ulasts.add(kmer_to_id[u[len(u) - k:]])

            line = reader.readline()
    return curr_id



def read_neikmer(neikmer, kmer_to_id, id_to_kmer, curr_id):
    '''
    :param neikmer:
    :return: all kmer here have >= startofblack id
    '''
    cnt = 0
    with open(neikmer, 'r') as reader:
        line = reader.readline()
        while (line != ""):
            kmer_to_append = line.rstrip()
            if not kmer_to_append in kmer_to_id:
                cnt += 1
                kmer_to_id[kmer_to_append] = curr_id
                id_to_kmer[curr_id] = kmer_to_append
                curr_id += 1
            line = reader.readline()
    return curr_id


def make_dbg(kmer_to_id, curr_id, dbg):
    suf_to_id = dict()
    pre_to_id = dict()
    for kmer, index in kmer_to_id.items():
        suf = kmer[1:k]
        pre = kmer[0:k - 1]
        if not suf in suf_to_id:
            suf_to_id[suf] = []
        suf_to_id[suf].append(index)

        if not pre in pre_to_id:
            pre_to_id[pre] = []
        pre_to_id[pre].append(index)

    for i in range(0, curr_id):
        kmer2 = id_to_kmer[i]
        suff = kmer2[1:k]
        prep = kmer2[0:k - 1]
        suflist = []
        prelist = []
        ins = False
        inp = False
        if suff in suf_to_id:
            ins = True
            suflist = suf_to_id[suff]
        if suff in pre_to_id:
            inp = True
            prelist = pre_to_id[suff]

        for j in suflist:
            for m in prelist:
                if not j in dbg:
                    dbg[j] = []
                dbg[j].append(m)
        if inp:
            del pre_to_id[suff]
        if ins:
            del suf_to_id[suff]

    # AT THIS POINT GRAPH IS Generated


def populate_in_out(id_to_kmer, indegree, outdegree):
    for key, value in id_to_kmer.items():
        if not key in dbg:
            dbg[key] = []
        indegree.append(0)
        outdegree.append(len(dbg[key]))
        # print(key, outdegree[len(outdegree)-1])

    for key, value in id_to_kmer.items():
        for j in dbg[key]:
            indegree[j] += 1

    # for key, value in id_to_kmer.items():
    #     print(key, "in=",indegree[key], ",out=",outdegree[key])


def read_copycount(copycount_file, copycount):
    with open(copycount_file, 'r') as reader:
        line = reader.readline()
        counter = 0
        while line != '':  # The EOF char is an empty string
            line = line.rstrip()
            copycount[counter] = int(line)
            counter += 1
            line = reader.readline()



def out_non_compacted_dbg_ve(id_to_kmer, dbg, VERTEX_FILE, EDGE_FILE):
    '''
    Print only vertex.txt and edges.txt
    :param id_to_kmer:
    :param dbg:
    :return:
    '''
    outvertices = open(VERTEX_FILE, 'w')
    for key, value in id_to_kmer.items():
        outvertices.write("{0} {1}\n".format(key, value))
    outvertices.close()

    outedges = open(EDGE_FILE, 'w')
    for key, value in dbg.items():
        for i in value:
            outedges.write("{0} {1}\n".format(key, i))
    outedges.close()

def out_non_compacted_dbg(id_to_kmer, dbg, startofblack, copycount, OLDGRAPH_FILE):
    '''
    GRAPH file print
    :param id_to_kmer:
    :param dbg:
    :param startofblack:
    :param copycount:
    :param OLDGRAPH_FILE:
    :return:
    '''
    out = open(OLDGRAPH_FILE, 'w')
    out.write("digraph d {\n")
    for key, value in id_to_kmer.items():
        # if key >= startofblack:
        #     out.write("{0} [label=\"{1}\", color=\"black\"]\n".format(key, id_to_kmer[key]))
        # elif(key in firsts and key in lasts):
        #     out.write("{0} [label=\"{1}\", color=\"purple\"]\n".format(key, id_to_kmer[key]))
        # elif(key in firsts and key not in lasts):
        #     out.write("{0} [label=\"{1}\", color=\"blue\"]\n".format(key, id_to_kmer[key]))
        # elif(key in lasts and key not in firsts):
        #     out.write("{0} [label=\"{1}\", color=\"red\"]\n".format(key, id_to_kmer[key]))
        # else:
        #     out.write("{0} [label=\"{1}\", color=\"green\"]\n".format(key, id_to_kmer[key]))

        if key >= startofblack:
            out.write("{0} [label=\"{1} ({2})\", color=\"black\"]\n".format(key, key, copycount[key]))
        elif (key in firsts and key in lasts):
            out.write("{0} [label=\"{1} ({2})\", color=\"purple\"]\n".format(key, key, copycount[key]))
        elif (key in firsts and key not in lasts):
            out.write("{0} [label=\"{1} ({2})\", color=\"blue\"]\n".format(key, key, copycount[key]))
        elif (key in lasts and key not in firsts):
            out.write("{0} [label=\"{1} ({2})\", color=\"red\"]\n".format(key, key, copycount[key]))
        else:
            out.write("{0} [label=\"{1} ({2})\", color=\"green\"]\n".format(key, key, copycount[key]))

    for key, value in dbg.items():
        for i in value:
            # out.write("{0} -> {1}[arrowhead=\"none\"]\n".format(key, i))
            out.write("{0} -> {1}[arrowhead=\"normal\"]\n".format(key, i))
    out.write("}\n")

    out.close()


def get_reverse_gr(dbg):
    revgr = dict()
    for key, value in dbg.items():
        for i in value:
            if i not in revgr:
                revgr[i] = []
            revgr[i].append(key)
    return revgr

def dfs(graph, start, parent_of, visited, group_id, group_assign, group_end, shape, group_len, unitigendlist):
    list1 = []
    list2 = []
    # parent_of = {}
    # visited = set()
    stack = [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            # list1.append(parent_of.get(vertex, None))
            list2.append(vertex)

            if (vertex in unitigendlist):
                group_assign[vertex] = group_id
                #print(vertex, id_to_kmer[vertex], end='\n')
                group_end[group_id] = vertex

                if (isstart_of_inu[vertex]):
                    shape[group_id] = "box"
                if (isend_of_inu[vertex]):
                    shape[group_id] = "diamond"

                group_len[group_id] = len(list2)
                return
            else:
                group_assign[vertex] = group_id

                if (isstart_of_inu[vertex]):
                    shape[group_id] = "box"
                if (isend_of_inu[vertex]):
                    shape[group_id] = "diamond"
                #print(vertex, id_to_kmer[vertex], end=' ')

            new_vertexes = set(graph[vertex]) - visited
            stack.extend(new_vertexes)
            # for i in new_vertexes:
            #     parent_of[i] = vertex
    # return parent_of


class CompactGraph:
    def __init__(self, compact_gr, group_end, group_assign, group_len, shape, unitigendlist):
        self.compact_gr = compact_gr
        self.group_end = group_end
        self.group_assign = group_assign
        self.group_len = group_len
        self.shape = shape
        self.unitigendlist = unitigendlist


def get_compact_graph(indegree, outdegree, genomefirst, genomelast, id_to_kmer, startofblack, firsts, lasts, ufirsts, ulasts, issampled, isstart_of_inu, isend_of_inu):
    # DO ADD GENOME FIRST AND LAST
    for key, value in genomefirst.items():
        # to make sure unitig first last are separately shown and not compacted
        if indegree[key] == 1:
            indegree[key] = 10
        if outdegree[key] == 1:
            outdegree[key] = 10

    for key in ufirsts:
        if indegree[key] == 1:
            indegree[key] = 10
        if outdegree[key] == 1:
            outdegree[key] = 10
    for key in ulasts:
        if indegree[key] == 1:
            indegree[key] = 10
        if outdegree[key] == 1:
            outdegree[key] = 10

    for key, value in genomelast.items():
        if indegree[key] == 1:
            indegree[key] = 10
        if outdegree[key] == 1:
            outdegree[key] = 10

    # findunitig start
    revgr = get_reverse_gr(dbg)
    unitigstartlist = []
    unitigendlist = []
    for key, value in id_to_kmer.items():
        if key >= startofblack:
            issampled[key] = False
            isstart_of_inu[key] = False
            isend_of_inu[key] = False
        elif (key in firsts and key in lasts):
            issampled[key] = True
            isstart_of_inu[key] = True
            isend_of_inu[key] = True
            # out.write("{0} [label=\"{1}\", color=\"purple\"]\n".format(key, key))
        elif (key in firsts and key not in lasts):
            issampled[key] = True
            isstart_of_inu[key] = True
            isend_of_inu[key] = False
            # out.write("{0} [label=\"{1}\", color=\"blue\"]\n".format(key, key))
        elif (key in lasts and key not in firsts):
            issampled[key] = True
            isstart_of_inu[key] = False
            isend_of_inu[key] = True
            # out.write("{0} [label=\"{1}\", color=\"red\"]\n".format(key, key))
        else:
            issampled[key] = True
            isstart_of_inu[key] = False
            isend_of_inu[key] = False
            # out.write("{0} [label=\"{1}\", color=\"green\"]\n".format(key, key))

        if indegree[key] == 1:
            parent = revgr[key][0]
            if (outdegree[parent] > 1):
                unitigstartlist.append(key)
        else:
            unitigstartlist.append(key)

        if outdegree[key] == 1:
            outparent = dbg[key][0]
            if (indegree[outparent] > 1):
                unitigendlist.append(key)
        else:
            unitigendlist.append(key)


    shape = {}
    group_len = {}
    group_assign = {}
    visited = set()
    parent_of = {}
    group_id = 0
    group_end = {}
    for key in unitigstartlist:
        if key not in visited:
            dfs(dbg, key, parent_of, visited, group_id, group_assign, group_end, shape, group_len, unitigendlist)
            group_id += 1
            # print(parent_of)
    # print(parent_of)
    compact_gr = {}
    for startvertex in unitigstartlist:
        grid = group_assign[startvertex]
        endvertex = group_end[grid]

        # for all in neigh of startvertex
        if startvertex not in revgr:
            revgr[startvertex] = []
        for inneigh in revgr[startvertex]:
            if group_assign[inneigh] not in compact_gr:
                compact_gr[group_assign[inneigh]] = set()
            compact_gr[group_assign[inneigh]].add(grid)

        if grid not in compact_gr:
            compact_gr[grid] = set()
        # for all out neigh of endvertex
        for outneigh in dbg[endvertex]:
            compact_gr[grid].add(group_assign[outneigh])

    cc_obj = CompactGraph(compact_gr, group_end, group_assign, group_len, shape, unitigendlist)
    return cc_obj


'''
print compact GRAPH
'''
def print_compact_graph(cc_obj, copycount):
    compact_gr = cc_obj.compact_gr
    group_end = cc_obj.group_end
    group_assign = cc_obj.group_assign
    group_len = cc_obj.group_len
    shape = cc_obj.shape

    compout = open("graph.gv", 'w')
    compout.write("digraph d {\n")
    compout.write("node [colorscheme=SVG] \n")

    iscolored = {}
    for key, value in compact_gr.items():
        # compout.write("{0} [label=\"{1}\", color=\"black\"]\n".format(key, key))
        endvertex = group_end[key]
        if endvertex not in vertex_to_genomeid:
            cc = 0
        else:
            arb_gen_id = next(iter(vertex_to_genomeid[endvertex]))
            if (len(vertex_to_genomeid[endvertex]) > 1):
                cc = 1
            else:
                cc = arb_gen_id % (len(colors)-2) + 2
            # cc=vertex_to_genomeid[endvertex] + 2

        shapp = "ellipse"
        if group_assign[endvertex] in shape:
            shapp = shape[group_assign[endvertex]]

        if endvertex in genomefirst:
            compout.write(
                "{0} [label=\"{1} ({2}){5}\", shape=\"{3}\", penwidth=5, style=\"dashed\", color=\"{4}\"]\n".format(key,
                                                                                                                    endvertex,
                                                                                                                    copycount[
                                                                                                                        endvertex],
                                                                                                                    shapp,
                                                                                                                    colors[
                                                                                                                        genomefirst[
                                                                                                                            endvertex]%(len(colors)-2) + 2],
                                                                                                                    group_len[
                                                                                                                        key]))
            iscolored[endvertex] = True
        if endvertex in genomelast:
            compout.write(
                "{0} [label=\"{1} ({2}){5}\", shape=\"{3}\", penwidth=5, style=\"dotted\", color=\"{4}\"]\n".format(key,
                                                                                                                    endvertex,
                                                                                                                    copycount[
                                                                                                                        endvertex],
                                                                                                                    shapp,
                                                                                                                    colors[
                                                                                                                        genomelast[
                                                                                                                            endvertex]%(len(colors)-2) + 2],
                                                                                                                    group_len[
                                                                                                                        key]))
            iscolored[endvertex] = True
        if endvertex in genomefirst and endvertex in genomelast and len(genomes[genomefirst[endvertex]]) != k:
            compout.write(
                "{0} [label=\"{1} ({2}){5}\", shape=\"{3}\", penwidth=5, style=\"diagonals\", color=\"{4}\"]\n".format(
                    key, endvertex, copycount[endvertex], shapp, colors[0], group_len[key]))
            iscolored[endvertex] = True
        if endvertex not in genomefirst and endvertex not in genomelast:
            compout.write(
                "{0} [label=\"{1} ({2}){5}\", shape=\"{3}\", penwidth=5, style=solid, color=\"{4}\"]\n".format(key,
                                                                                                               endvertex,
                                                                                                               copycount[
                                                                                                                   endvertex],
                                                                                                               shapp,
                                                                                                               colors[
                                                                                                                   cc],
                                                                                                               group_len[
                                                                                                                   key]))
            iscolored[endvertex] = True

    for key, value in compact_gr.items():
        for i in value:
            # out.write("{0} -> {1}[arrowhead=\"none\"]\n".format(key, i))
            compout.write("{0} -> {1}[arrowhead=\"normal\"]\n".format(key, i))
    compout.write("}\n")
    compout.close()


def print_compact_graph_misjoin(cc_obj, copycount, vertices_to_print, GRAPH_MISJOIN_FILE):
    compact_gr = cc_obj.compact_gr
    group_end = cc_obj.group_end
    group_assign = cc_obj.group_assign
    group_len = cc_obj.group_len
    shape = cc_obj.shape

    compout = open(GRAPH_MISJOIN_FILE, 'w')
    compout.write("digraph d {\n")
    compout.write("node [colorscheme=SVG] \n")

    iscolored = {}
    for key, value in compact_gr.items():
        if key not in vertices_to_print:
            continue
        # compout.write("{0} [label=\"{1}\", color=\"black\"]\n".format(key, key))
        endvertex = group_end[key]
        if endvertex not in vertex_to_genomeid:
            cc = 0
        else:
            arb_gen_id = next(iter(vertex_to_genomeid[endvertex]))
            if (len(vertex_to_genomeid[endvertex]) > 1):
                cc = 1
            else:
                cc = (arb_gen_id % (len(colors)-2)) + 2
            # cc=vertex_to_genomeid[endvertex] + 2

        shapp = "ellipse"
        if group_assign[endvertex] in shape:
            shapp = shape[group_assign[endvertex]]

        if endvertex in genomefirst:
            compout.write(
                "{0} [label=\"{1} ({2}){5}\", shape=\"{3}\", penwidth=5, style=\"dashed\", color=\"{4}\"]\n".format(key,
                                                                                                                    endvertex,
                                                                                                                    copycount[
                                                                                                                        endvertex],
                                                                                                                    shapp,
                                                                                                                    colors[
                                                                                                                        genomefirst[
                                                                                                                            endvertex] %(len(colors)-2) + 2],
                                                                                                                    group_len[
                                                                                                                        key]))
            iscolored[endvertex] = True
        if endvertex in genomelast:
            compout.write(
                "{0} [label=\"{1} ({2}){5}\", shape=\"{3}\", penwidth=5, style=\"dotted\", color=\"{4}\"]\n".format(key,
                                                                                                                    endvertex,
                                                                                                                    copycount[
                                                                                                                        endvertex],
                                                                                                                    shapp,
                                                                                                                    colors[
                                                                                                                        genomelast[
                                                                                                                            endvertex] %(len(colors)-2) + 2],
                                                                                                                    group_len[
                                                                                                                        key]))
            iscolored[endvertex] = True
        if endvertex in genomefirst and endvertex in genomelast and len(genomes[genomefirst[endvertex]]) != k:
            compout.write(
                "{0} [label=\"{1} ({2}){5}\", shape=\"{3}\", penwidth=5, style=\"diagonals\", color=\"{4}\"]\n".format(
                    key, endvertex, copycount[endvertex], shapp, colors[0], group_len[key]))
            iscolored[endvertex] = True
        if endvertex not in genomefirst and endvertex not in genomelast:
            compout.write(
                "{0} [label=\"{1} ({2}){5}\", shape=\"{3}\", penwidth=5, style=solid, color=\"{4}\"]\n".format(key,
                                                                                                               endvertex,
                                                                                                               copycount[
                                                                                                                   endvertex],
                                                                                                               shapp,
                                                                                                               colors[
                                                                                                                   cc],
                                                                                                               group_len[
                                                                                                                   key]))
            iscolored[endvertex] = True

    for key, value in compact_gr.items():
        if key not in vertices_to_print:
            continue
        for i in value:
            # out.write("{0} -> {1}[arrowhead=\"none\"]\n".format(key, i))
            compout.write("{0} -> {1}[arrowhead=\"normal\"]\n".format(key, i))
    compout.write("}\n")
    compout.close()

def write_html(GTOUCH_FILE, HTML_FILE):
    out = open(HTML_FILE, "w")
    out.write("<html>")
    out.write("<body>")

    with open(GTOUCH_FILE, 'r') as reader:
        line = reader.readline()
        out.write("<p>")
        oldinvalid = ""
        flag = False
        if line == '':
            out.write("</p>")
        while line != '':  # The EOF char is an empty string
            line = line.rstrip()
            linesplit = line.split()
            invaliduni = linesplit[0].rstrip()
            if (invaliduni != oldinvalid):
                flag = True
            if (flag):
                out.write("</p>")
                out.write("<p>")
                bb = ""
                rr = ""
                pp = ""
                kk = ""
                if (len(invaliduni) >= 2 * k):
                    bb = invaliduni[0:k]
                    kk = invaliduni[k:len(invaliduni) - k]
                    rr = invaliduni[len(invaliduni) - k:]
                else:
                    bb = invaliduni[0:len(invaliduni) - k]
                    pp = invaliduni[len(invaliduni) - k:k]
                    rr = invaliduni[k:len(invaliduni)]
                # out.write("<u>"+invaliduni+"</u><br>")
                out.write(
                    "<u><span style='color:blue;'>{0}</span><span style='color:black;'>{1}</span><span style='color:purple;'>{2}</span><span style='color:red;'>{3}</span></u><br>".format(
                        bb, kk, pp, rr))
                flag = False
            oldinvalid = invaliduni

            g = linesplit[1].rstrip()
            firstkmer = invaliduni[0:k]
            lastkmer = invaliduni[len(invaliduni) - k:]

            line = reader.readline()

            l_start = g.find(lastkmer)
            f_start = g.find(firstkmer)

            assert (
                        l_start == -1 or f_start == -1 or l_start <= f_start), "u_m should appear before u_0 in an invalid unitig!"

            bluetext = ""
            redtext = ""
            purpletext = ""
            p1 = ""
            p2 = ""
            p3 = ""
            # f_start+k-1 >= l_start
            if (f_start != -1 and l_start != -1 and f_start <= l_start + k - 1):
                p1 = g[0:l_start]
                redtext = g[l_start:f_start]
                purpletext = g[f_start:l_start + k]
                bluetext = g[l_start + k:f_start + k]
                p3 = g[f_start + k:]
            elif (f_start != -1 and l_start != -1 and f_start > l_start + k - 1):
                p1 = g[0:l_start]
                redtext = g[l_start:l_start + k]
                p2 = g[l_start + k:f_start]
                bluetext = g[f_start:f_start + k]
                p3 = g[f_start + k:]
            elif (f_start != -1 and l_start == -1):
                p1 = g[0:f_start]
                bluetext = firstkmer
                p3 = g[f_start + k:]
            elif (f_start == -1 and l_start != -1):
                p1 = g[0:l_start]
                redtext = lastkmer
                p3 = g[l_start + k:]
            else:
                p1 = g
            out.write(
                "<span style='color:{6};'>{0}</span><span style='color:red;'>{1}</span><span style='color:purple;'>{2}</span><span style='color:{6};'>{3}</span><span style='color:blue;'>{4}</span><span style='color:{6};'>{5}</span><br>".format(
                    p1, redtext, purpletext, p2, bluetext, p3, "black"))
            # colors[vertex_to_genomeid[allunitigkmers[lastkmer]]]
    # 0+1+        2+      3   4       5
    # p1+redtext+purpletext+p2+bluetext+p3

    out.write("</body>")
    out.write("</html>")

    out.close()


class GraphSCC:
    # init function to declare class variables
    def __init__(self, compact_gr):
        self.V = len(compact_gr)
        self.adj = [[] for i in range(self.V)]
        for key, value in compact_gr.items():  # add all edges (unidirected)
            for i in value:
                self.adj[key].append(i)
                self.adj[i].append(key)

    def DFSUtil(self, temp, v, visited):
        # Mark the current vertex as visited
        visited[v] = True

        # Store the vertex to list
        temp.append(v)

        # Repeat for all vertices adjacent
        # to this vertex v
        for i in self.adj[v]:
            if visited[i] == False:
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp

    # Method to retrieve connected components
    # in an undirected graph
    def connectedComponents(self, list_of_vertices):
        visited = []
        cc = []
        for i in range(self.V):
            visited.append(False)
        for v in list_of_vertices:
            if visited[v] == False:
                temp = []
                for all in self.DFSUtil(temp, v, visited):
                    cc.append(all)
        return cc

def get_vertices_in_misjoin(dbg,unitigs,cc_obj,kmer_to_id,issampled):
    revdbg=get_reverse_gr(dbg)
    list_of_unitig_firsts=[]
    for u in unitigs:
        ui = -1
        uj = -1

        u0 = u[0:k]
        un = u[len(u) - k:]

        # if not (copycount[kmer_to_id[u0]]==1 and copycount[kmer_to_id[u0]]==1):
        #     list_of_unitig_firsts.append(cc_obj.group_assign[kmer_to_id[u0]])

        iprimes = dict()  # start of some genome
        jprimes = dict()  # end of some genome
        for i in range(len(u) - k + 1):
            kmer=u[i:i + k]
            kmerid=kmer_to_id[kmer]
            ind=len(revdbg[kmer_to_id[kmer]])
            outd=len(dbg[kmer_to_id[kmer]])

            if (ind == 2 and outd == 1):
                if (copycount[kmerid] == 2):
                    for inneigh in revdbg[kmerid]:
                        if (not issampled[inneigh]):
                            if (copycount[inneigh] == 1):
                                ui = i

            if (ind == 1 and outd == 1):
                if (copycount[kmer_to_id[kmer]]!=2):
                    ui=-1

            # if(ind==2):
            #     if (copycount[kmerid]!=2):
            #         ui=-1
            #         continue

            if (outd==2 and ind == 1):
                if (copycount[kmerid]==2):
                    for outneigh in dbg[kmer_to_id[kmer]]:
                        if (not issampled[outneigh]):
                            if (copycount[outneigh] == 1):
                                if(ui!=-1):
                                    uj = i

            # if(ind==2 and outd==2):
            #     print("Error: ")
                #exit(1)

            if(i>=ui and ui!=-1):
                if(kmerid in genomefirst):
                    if i not in iprimes:
                        iprimes[i]=[]
                    iprimes[i].append(genomefirst[kmerid])
                if (kmerid in genomelast):
                    if i not in jprimes:
                        jprimes[i]=[]
                    jprimes[i].append(genomelast[kmerid])

            check=False
            if(ui!=-1 and uj!=-1):
                for iprime, gi in iprimes.items():
                    for jprime, gj in jprimes.items():
                        # gi=iprimes[iprime]
                        # gj=jprimes[jprime]
                        check1 = (iprime-1 <= jprime)
                        check2 = (gi!=gj)
                        check = (check1 and check2)
                        if(check):
                            break
                    if (check):
                        break
                if(check):
                    assert(jprime<=uj)
                    assert (iprime <= uj)
                    list_of_unitig_firsts.append(cc_obj.group_assign[kmer_to_id[u[ui:ui + k]]])
                    break
                else:
                    ui=-1
                    uj=-1



    print("Total true invalid:", len(unitigs))
    print("True invalid unitig containing misjoin", len(list_of_unitig_firsts))
    print("True invalid unitig containing NO misjoin", len(unitigs)-len(list_of_unitig_firsts))
    scc=GraphSCC(cc_obj.compact_gr)
    vertices_from_misjoin=set(scc.connectedComponents(list_of_unitig_firsts))
    # vertices_to_print=set()
    # for i in range(len(dbg)):
    #     if i not in vertices_from_misjoin:
    #         vertices_to_print.add(i)
    return vertices_from_misjoin



def get_vertices_in_misjoin_plus(dbg,unitigs,cc_obj,kmer_to_id,issampled):
    revdbg=get_reverse_gr(dbg)
    list_of_unitig_firsts=[]

    for u in unitigs:

        setkmerid = -1
        ui = -1
        uj = -1

        u0 = u[0:k]
        un = u[len(u) - k:]

        # if not (copycount[kmer_to_id[u0]]==1 and copycount[kmer_to_id[u0]]==1):
        #     list_of_unitig_firsts.append(cc_obj.group_assign[kmer_to_id[u0]])

        iprimes = dict()  # start of some genome
        jprimes = dict()  # end of some genome
        for i in range(len(u) - k + 1):
            kmer=u[i:i + k]
            kmerid=kmer_to_id[kmer]
            ind=len(revdbg[kmer_to_id[kmer]])
            outd=len(dbg[kmer_to_id[kmer]])

            #if (ind == 2 and outd == 1):#
            if (outd == 1):  #
                if (copycount[kmerid] >= 2): #2
                    for inneigh in revdbg[kmerid]:
                        if (not issampled[inneigh]):
                            #if (copycount[inneigh] == copycount[kmerid]-1): #1
                            ui = i
                            setkmerid=copycount[kmerid]


            if (ind == 1 and outd == 1):
                if (copycount[kmer_to_id[kmer]]!=setkmerid): #2
                    ui=-1

            # if(ind==2):
            #     if (copycount[kmerid]!=2):
            #         ui=-1
            #         continue

            if (outd>=2 and ind == 1):
                if (copycount[kmerid]==setkmerid): #2
                    for outneigh in dbg[kmer_to_id[kmer]]:
                        if (not issampled[outneigh]):
                            #if (copycount[outneigh] == setkmerid-1): #1
                            if(ui!=-1):
                                uj = i

            # if(ind==2 and outd==2):
            #     print("Error: ")
                #exit(1)

            if(i>=ui and ui!=-1):
                if(kmerid in genomefirst):
                    if i not in iprimes:
                        iprimes[i]=[]
                    iprimes[i].append(genomefirst[kmerid])
                if (kmerid in genomelast):
                    if i not in jprimes:
                        jprimes[i]=[]
                    jprimes[i].append(genomelast[kmerid])

            check=False
            if(ui!=-1 and uj!=-1):
                for iprime, gi in iprimes.items():
                    for jprime, gj in jprimes.items():
                        # gi=iprimes[iprime]
                        # gj=jprimes[jprime]
                        check1 = (iprime-1 <= jprime)
                        check2 = (gi!=gj)
                        check = (check1 and check2)
                        if(check):
                            break
                    if (check):
                        break
                if(check):
                    assert(jprime<=uj)
                    assert (iprime <= uj)
                    list_of_unitig_firsts.append(cc_obj.group_assign[kmer_to_id[u[ui:ui + k]]])
                    break
                else:
                    ui=-1
                    uj=-1



    print("Total true invalid:", len(unitigs))
    print("True invalid unitig containing misjoin", len(list_of_unitig_firsts))
    print("True invalid unitig containing NO misjoin", len(unitigs)-len(list_of_unitig_firsts))
    scc=GraphSCC(cc_obj.compact_gr)
    vertices_from_misjoin=set(scc.connectedComponents(list_of_unitig_firsts))
    # vertices_to_print=set()
    # for i in range(len(dbg)):
    #     if i not in vertices_from_misjoin:
    #         vertices_to_print.add(i)
    return vertices_from_misjoin





def write_copycount_jf_to_file(vertex_f, copycount_f):
    #sset_to_fa.sh ref.sset > ref.fa
    #jellyfish count -m $K -s 100M -t 10 ref.fa
    os.system("cat {0} | cut -f2 -d\" \" | jellyfish query -i mer_counts.jf  | cut -f2 -d \" \" > {1}".format(vertex_f, copycount_f))


def write_one_hop_from_iu(unitigs, k, outfile):
    invalidnei = open(outfile, 'w')
    for invaliduni in unitigs:
        for i in range(len(invaliduni) - k + 1):
            toappend = invaliduni[i:i + k]
            invalidnei.write("A" + toappend[0:k - 1] + "\n")
            invalidnei.write("C" + toappend[0:k - 1] + "\n")
            invalidnei.write("G" + toappend[0:k - 1] + "\n")
            invalidnei.write("T" + toappend[0:k - 1] + "\n")
            invalidnei.write(toappend[1:k] + "A\n")
            invalidnei.write(toappend[1:k] + "C\n")
            invalidnei.write(toappend[1:k] + "G\n")
            invalidnei.write(toappend[1:k] + "T\n")
    invalidnei.close()


def write_neikmer_jf_to_file(NEIKMER_FILE):
    write_one_hop_from_iu(unitigs, k, "invalidnei1.txt")
    returned_value = os.system("cat invalidnei1.txt | jellyfish query -i mer_counts.jf > neikmer1.tmp")
    # os.system("cat invalidnei1.txt | jellyfish query -i mer_counts.jf > neikmer1.tmp")
    os.system("paste -d \" \" invalidnei1.txt neikmer1.tmp | awk '$2 != \"0\"' | cut -f 1 -d \" \" > {0}".format(NEIKMER_FILE))
    os.system("rm neikmer1.tmp invalidnei1.txt")

def dot_vis():
    #python dbgofinvalid.py gtouch 21;
    os.system("neato -Tsvg graph_misjoin.gv -o dbg_mis.svg")
    os.system("neato -Tsvg graph_no_misjoin.gv -o dbg_nomis.svg")
    os.system("neato -Tsvg graph.gv -o dbg.svg")
    os.system("neato -Tsvg oldgraph.gv -o olddbg.svg")
#    os.system("neato -Tsvg graph.gv -o dbg.svg; neato -Tsvg oldgraph.gv -o olddbg.svg; neato -Tsvg graph_misjoin.gv -o dbg_mis.svg")
#   os.system("neato -Tsvg graph.gv -o dbg.svg; neato -Tsvg oldgraph.gv -o olddbg.svg")



kmer_to_id = dict()  # key=kmer value=kmer-id
id_to_kmer = dict()  # key=kmer-id value=kmer
dbg = dict()

indegree = []
outdegree = []

firsts = []  # kmer-id
lasts = []  # kmer-id

ufirsts = set()  # kmer-id
ulasts =  set() # kmer-id

genomefirst = {}  # key=kmer, value=genome-id
genomelast = {}  # key=kmer, value=genome-id

genome_to_gid = {}  ##key=genomestring, value=genomeid
genomes = []  # strings genome
vertex_to_genomeid = {}  ###key=kmer, value=genomeids[]

issampled = {}
isstart_of_inu = {}
isend_of_inu = {}

copycount = {}

unitigs=[]

NEIKMER_FILE="neikmer.txt"
VERTEX_FILE="vertices.txt"
EDGE_FILE="edges.txt"
COPYCOUNT_FILE="copycount.txt"
OLDGRAPH_FILE="oldgraph.gv"
HTML_FILE="out.html"

#curr_id = read_gtouch(GTOUCH_FILE, kmer_to_id, id_to_kmer, genome_to_gid, genomes, genomefirst, genomelast, firsts, lasts, ufirsts, ulasts, unitigs)
curr_id = read_gtouch_idversion(INVALID_UNITIG_FILE, GENOME_FILE, TOUCHED_GENOME_ID_FILE, kmer_to_id, id_to_kmer, genome_to_gid, genomes, genomefirst, genomelast, firsts, lasts, ufirsts, ulasts, unitigs)
startofblack = curr_id
write_neikmer_jf_to_file(NEIKMER_FILE)
curr_id = read_neikmer(NEIKMER_FILE, kmer_to_id, id_to_kmer, curr_id)
make_dbg(kmer_to_id, curr_id, dbg)
out_non_compacted_dbg_ve(id_to_kmer, dbg, VERTEX_FILE, EDGE_FILE)

populate_in_out(id_to_kmer, indegree, outdegree)

write_copycount_jf_to_file(VERTEX_FILE, COPYCOUNT_FILE)
read_copycount(COPYCOUNT_FILE, copycount)
out_non_compacted_dbg(id_to_kmer, dbg, startofblack, copycount, OLDGRAPH_FILE)

## NOW COMPACTED ONE
cc_obj=get_compact_graph(indegree, outdegree, genomefirst, genomelast, id_to_kmer, startofblack, firsts, lasts, ufirsts, ulasts, issampled, isstart_of_inu, isend_of_inu)
print_compact_graph(cc_obj, copycount)
#write_html(GTOUCH_FILE, HTML_FILE)

vertices_from_misjoin=get_vertices_in_misjoin_plus(dbg,unitigs,cc_obj,kmer_to_id,issampled)
vertices_not_in_misjoin=set()
for i in range(len(dbg)):
    if i not in vertices_from_misjoin:
        vertices_not_in_misjoin.add(i)

print_compact_graph_misjoin(cc_obj, copycount, vertices_from_misjoin, "graph_misjoin.gv")
print_compact_graph_misjoin(cc_obj, copycount, vertices_not_in_misjoin, "graph_no_misjoin.gv")
dot_vis()



# repst = open("repst", "w")
# repst.write(str(startofblack))
# repst.close()






