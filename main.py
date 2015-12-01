__author__ = 'JDima'
import random
from numpy import random
import math
from functools import partial

reference = [23.077, 2.770, 1.506, 27.417, 3.242, 7.507, 6.510, 17.395, 1.282, 10.298, 18.440, 16.291, 1.297, 3.111, 12.021, 5.108, 5.622, 17.695, 22.296,
             14.167,	6.554,	9.762,	0.078,	12.606,	0.102,	2.339,	2.307,	1.853,	57.148,	1.450,	4.940,	0.388,	2.295,	8.115,	3.957,	2.387,	7.966,
             60.127,	0.921,	5.609,	8.514,	3.502,	0.323,	2.224,	7.729,	60.689,	24.140,	5.976,	0.196,	6.220,	4.057,	59.424,	8.143,	8.374,	1.795,	5.363,
             3.161	,0.123,	7.657,	0.637,	3.558,	56.790,	11.925,	16.404,	2.962,	3.932,	8.720,	10.715,	22.930,	2.963,	0.690,	3.448,	31.847,	20.789,	4.341,
             16.539	,63.213,	1.952,	0.034,	2.332,	0.189,	3.570,	22.823,	0.256,	7.688,	4.189,	0.661,	0.067,	2.157,	69.629,	15.588,	17.693,	36.930,	7.544,
             52.872	,9.977	,5.122,	13.237,	64.086,	29.131,	12.586,	12.287,	31.099,	5.166,	6.078,	6.060,	4.082,	4.989,	1.027,	6.845,	1.305,	8.560,	11.630,
             15.450	,8.022	,6.255,	3.768,	2.092,	0.220,	2.046,	9.987,	67.034,	0.772,	10.791,	0.777,	6.782,	4.640,	2.031,	16.929,	14.122,	3.754,	13.863,
             12.003	,0.034	,2.324,	0.188,	3.563,	6.468,	31.100,	12.313,	12.304,	31.211,	6.648,	59.863,	23.394,	0.392,	3.449,	3.181,	36.432,	2.695]

class Map:
    def __init__(self, cuts, ids):
        self.cuts = cuts
        self.ids = ids
        self.num_ids = {id : i for i, id in enumerate(ids)}

    def __str__(self):
        return "Frags: {}\nids: {}".format(list(map(round, self.cuts)), self.ids)




class MapManager:

    def __init__(self, n, sigma, pc):
        self.sigma = sigma
        self.pc = pc
        self.n = n
        self.map1, self.map2 = self.create_maps()

    def get_frag(self, icut):
        return reference[icut] + random.normal(0, self.sigma, 1)[0]

    def create_map(self):
        ids, cuts = [], [0]
        #which was added by pc
        added = 0
        for icut in range(self.n):
            cur_icut = icut
            icut += added
            #always add cut
            cur_cut = self.get_frag(icut)
            #save last added cut, it need when we repeat by pc

            for r in random.binomial(1, 1 - self.pc, len(reference)).tolist():
                if not r:
                    break
                added += 1
                icut += 1
                #save repeated id help to next calculations
                cur_cut += self.get_frag(icut)

            cuts.append(cuts[-1] + cur_cut)
            ids.append(icut)

        ids.append(icut + added + 1)
        return Map(cuts, ids)

    def create_maps(self):
        return self.create_map(), self.create_map()


class Experiment:

    def __init__(self, n, sigma, pc, gap):
        self.sigma = sigma
        self.gap = gap
        self.pc = pc
        self.n = n

    def get_pairs_in_frags(self, map1, map2):
        id2 = 0
        pairs = []
        for id1 in map1.ids[self.gap :]:
            # print(map2.ids[id2], id1)
            pairs.append((map2.ids[id2], id1))
            id2 += 1
        return pairs

    def get_pairs_in_kb(self, map1, map2, pairs):
        return [(map2.cuts[map2.num_ids[p2]], map1.cuts[map1.num_ids[p1]]) for p2, p1 in pairs]

    def calc_distances(self, pairs):
        return [math.fabs(x - y) / math.sqrt(2) for x, y in pairs]

    def get_in_border(self, sigma_dist, border_dist):
        reference = [1, 2, 3, 4, 5]
        for dist_sigma, dist_point in zip(sigma_dist, border_dist):
            ref_sum = 0
            count_frags_ref = 0
            for r in reference:
                ref_sum += r
                if dist_sigma >= ref_sum:
                    count_frags_ref += 1
                    continue
                else:
                    break
            if dist_point < math.sqrt(count_frags_ref) * 3 * self.sigma:
                return True
        return False


    def sigma_dist(self, pairs):
        # Проекцируем на ось, находим расстояние и отображаем ось х или у => ((x+y) / 2, (x+y) / 2) => d = (x+y)/sqrt(2) => отображение на sqrt(2) / 2 => имеем (x + y)/2
        return [(x + y) / 2 for x, y in pairs]

    def test_border(self, map1, map2):
        #Build pairs of border points
        pairs_in_frags = self.get_pairs_in_frags(map1, map2)

        # print(self.get_pairs_in_kb(map1, map2, left_pairs))
        #Recalc pairs in frags to pairs in kb
        pairs_kb = self.get_pairs_in_kb(map1, map2, pairs_in_frags)

        #calc array of distance from (0,0) to line in kb
        sigma_dist = self.sigma_dist(pairs_kb)

        #array of distances from obrder points to ideal align
        border_dist = self.calc_distances(pairs_kb)
        return self.get_in_border(sigma_dist, border_dist)

    def run(self):
        mp = MapManager(self.n, self.sigma, self.pc)
        return self.test_border(mp.map1, mp.map2) and \
               self.test_border(mp.map2, mp.map1)


# print(exp.run())

# f = open('analyze_gap.txt', 'w')
#
# N = 500
# for isigma in range(1,10):
#     sigma = isigma * 0.2
#     for ipc in range(1, 11):
#         pc = 0.5 + ipc * 0.05
#         for gap in range(0, 7):
#             exp = Experiment(25, sigma, pc, gap)
#             errors = 0
#             for _ in range(N):
#                 errors += exp.run()
#             print("{} {} {} {}".format(sigma, pc, gap, errors / N))
#             print("{} {} {} {}".format(sigma, pc, gap, errors / N), file=f)
# f.close()


exp = Experiment(25, 1.8, 0.6, 2)
map1 = Map([0, 1, 2, 5], [0, 1, 2, 5])
map2 = Map([0, 1, 3, 4, 5], [0, 1, 3, 4, 5])
print(exp.test_border(map1, map2))
