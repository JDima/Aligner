__author__ = 'JDima'
import math
import abc
from numpy import random
from multiprocessing.dummy import Pool as ThreadPool

reference_frags = [23.077, 2.770, 1.506, 27.417, 3.242, 7.507, 6.510, 17.395, 1.282, 10.298, 18.440, 16.291, 1.297,
                   3.111, 12.021, 5.108, 5.622, 17.695, 22.296,
                   14.167, 6.554, 9.762, 0.078, 12.606, 0.102, 2.339, 2.307, 1.853, 57.148, 1.450, 4.940, 0.388, 2.295,
                   8.115, 3.957, 2.387, 7.966,
                   60.127, 0.921, 5.609, 8.514, 3.502, 0.323, 2.224, 7.729, 60.689, 24.140, 5.976, 0.196, 6.220, 4.057,
                   59.424, 8.143, 8.374, 1.795, 5.363,
                   3.161, 0.123, 7.657, 0.637, 3.558, 56.790, 11.925, 16.404, 2.962, 3.932, 8.720, 10.715, 22.930,
                   2.963, 0.690, 3.448, 31.847, 20.789, 4.341,
                   16.539, 63.213, 1.952, 0.034, 2.332, 0.189, 3.570, 22.823, 0.256, 7.688, 4.189, 0.661, 0.067, 2.157,
                   69.629, 15.588, 17.693, 36.930, 7.544,
                   52.872, 9.977, 5.122, 13.237, 64.086, 29.131, 12.586, 12.287, 31.099, 5.166, 6.078, 6.060, 4.082,
                   4.989, 1.027, 6.845, 1.305, 8.560, 11.630,
                   15.450, 8.022, 6.255, 3.768, 2.092, 0.220, 2.046, 9.987, 67.034, 0.772, 10.791, 0.777, 6.782, 4.640,
                   2.031, 16.929, 14.122, 3.754, 13.863,
                   12.003, 0.034, 2.324, 0.188, 3.563, 6.468, 31.100, 12.313, 12.304, 31.211, 6.648, 59.863, 23.394,
                   0.392, 3.449, 3.181, 36.432, 2.695]

reference_dist = [0, 23.077, 25.847, 27.353, 54.77, 58.012, 65.519, 72.029, 89.424, 90.706, 101.004, 119.444, 135.735,
                  137.032, 140.143, 152.164, 157.272, 162.894, 180.589, 202.885, 217.052, 223.606, 233.368, 233.446,
                  246.052, 246.154, 248.493, 250.8, 252.653, 309.801, 311.251, 316.191, 316.579, 318.874, 326.989,
                  330.946, 333.333, 341.299, 401.426, 402.347, 407.956, 416.47, 419.972, 420.295, 422.519, 430.248,
                  490.937, 515.077, 521.053, 521.249, 527.469, 531.526, 590.95, 599.093, 607.467, 609.262, 614.625,
                  617.786, 617.909, 625.566, 626.203, 629.761, 686.551, 698.476, 714.88, 717.842, 721.774, 730.494,
                  741.209, 764.139, 767.102, 767.792, 771.24, 803.087, 823.876, 828.217, 844.756, 907.969, 909.921,
                  909.955, 912.287, 912.476, 916.046, 938.869, 939.125, 946.813, 951.002, 951.663, 951.73, 953.887,
                  1023.516, 1039.104, 1056.797, 1093.727, 1101.271, 1154.143, 1164.12, 1169.242, 1182.479, 1246.565,
                  1275.696, 1288.282, 1300.569, 1331.668, 1336.834, 1342.912, 1348.972, 1353.054, 1358.043, 1359.07,
                  1365.915, 1367.22, 1375.78, 1387.41, 1402.86, 1410.882, 1417.137, 1420.905, 1422.997, 1423.217,
                  1425.263, 1435.25, 1502.284, 1503.056, 1513.847, 1514.624, 1521.406, 1526.046, 1528.077, 1545.006,
                  1559.128, 1562.882, 1576.745, 1588.748, 1588.782, 1591.106, 1591.294, 1594.857, 1601.325, 1632.425,
                  1644.738, 1657.042, 1688.253, 1694.901, 1754.764, 1778.158, 1778.55, 1781.999, 1785.18, 1821.612,
                  1824.307]

# class BinarySearch:
#
#     def __init__(self, task):
#         self.task = task
#
#     def binarySearch(self, left, right):
#         if right < left:
#             return self.NOT_FOUND
#
        # mid = (left + right) / 2
        # task_result = self.task(mid)
        # if task_result > 0:
        #     return self.binarySearch(arr, searchValue, mid + 1, right)
        # elif task_result < 0:
        #     return self.binarySearch(arr, searchValue, left, mid - 1)
        # else:
        #     return mid
#
#     def search(self, left, right, ):
#         self.left = left
#         self.right = right
#         return self.binarySearch(left, right)
#




class Map:
    def __init__(self, cuts, ids):
        self.cuts = cuts
        self.ids = ids
        self.id_nums = {id: i for i, id in enumerate(ids)}
        self.num_ids = {i: id for i, id in enumerate(ids)}

    def __str__(self):
        return "Frags: {}\nids: {}".format(list(map(round, self.cuts)), self.ids)


class MapManager:
    def __init__(self, n, sigma, pc):
        self.sigma = sigma
        self.pc = pc
        self.n = n
        self.map1, self.map2 = self.create_maps()

    def get_frag(self, icut):
        return reference_frags[icut] + random.normal(0, self.sigma, 1)[0]

    def create_map(self):
        ids, cuts = [0], [0]
        # which was added by pc
        added = 0
        for icut in range(self.n):
            icut += added
            if icut >= self.n :
                break
            # always add cut
            cur_cut = reference_frags[icut]
            # save last added cut, it need when we repeat by pc

            for r in random.binomial(1, 1 - self.pc, len(reference_frags)).tolist():
                if not r:
                    break
                added += 1
                icut += 1
                # save repeated id help to next calculations
                cur_cut += reference_frags[icut]
                if icut >= self.n:
                    break

            cuts.append(cuts[-1] + cur_cut + random.normal(0, self.sigma, 1)[0])
            ids.append(icut + 1)

        # ids.append(icut + added + 1)
        return Map(cuts, ids)

    def create_maps(self):
        return self.create_map(), self.create_map()



class Experiment:
    def __init__(self, map_size, sigma, pc, gap):
        self.sigma = sigma
        self.gap = gap
        self.pc = pc
        self.map_size = map_size

    @abc.abstractmethod
    def run(self):
        pass

    @staticmethod
    def get_filename():
        pass

class FragmentExperiment(Experiment):

    def get_pairs_in_frags(self, map1, map2):
        num_id2 = 0
        pairs = {}
        for id1 in map1.ids[self.gap:]:
            if num_id2 not in map2.num_ids:
                break
            pairs[map2.num_ids[num_id2]] = id1
            num_id2 += 1
        return pairs

    def out_border(self, map1, map2):
        alignment_ids = sorted(list(set(map1.ids) & set(map2.ids)))
        pairs_cuts = self.get_pairs_in_frags(map1, map2)

        for align_id in alignment_ids:
            if pairs_cuts.get(align_id, 10**10) < align_id:
                return True
        return False

    def run(self):
        mp = MapManager(self.map_size, self.sigma, self.pc)
        return self.out_border(mp.map1, mp.map2) and \
               self.out_border(mp.map2, mp.map1)

    @staticmethod
    def get_filename():
        return "gap_in_frag.txt"


class FragmentExperimentSize(Experiment):

    def get_pairs_in_frags(self, map1, map2):
        num_id2 = 0
        pairs = {}
        for id1 in map1.ids[self.gap:]:
            if num_id2 not in map2.num_ids:
                break
            pairs[map2.num_ids[num_id2]] = id1
            num_id2 += 1
        return pairs


    def invert(self, dictionary):
        inverted_dict = {value: 0 for value in dictionary.values()}
        for key in dictionary:
            inverted_dict[dictionary[key]] = key
        return inverted_dict

    def count_in_region(self, map1, map2):
        pairs_cuts_down = self.get_pairs_in_frags(map1, map2)
        pairs_cuts_up = self.invert(self.get_pairs_in_frags(map2, map1))

        count_in_region = 0
        for cut in map2.ids:
            up = pairs_cuts_up.get(cut, 0)
            down = pairs_cuts_down.get(cut, map1.ids[-1])
            count_in_region += map1.id_nums[down] - map1.id_nums[up] + 1
            if down == map1.ids[-1] == up:
                break

        return count_in_region / (len(map1.cuts) * len(map2.cuts))


    def run(self):
        mp = MapManager(self.map_size, self.sigma, self.pc)
        return self.count_in_region(mp.map2, mp.map1)

    @staticmethod
    def get_filename():
        return "gap_in_frag_size.txt"


class KilobaseExperimentSize(Experiment):
    def __init__(self, map_size, sigma, pc, gap):
        Experiment.__init__(self, map_size, sigma, pc, gap)
        # self.gap *= self.sigma

    def count_in_region(self, map1, map2):
        l = (len(map1.cuts) + len(map2.cuts)) / 2
        count_in_region = 0
        for cut1 in map1.cuts:
            for cut2 in map2.cuts:
                if -self.gap * math.sqrt(l) <= cut2 - cut1 <= self.gap * math.sqrt(l):
                    count_in_region += 1
        return count_in_region / (len(map1.cuts) * len(map2.cuts))


    def run(self):
        mp = MapManager(self.map_size, self.sigma, self.pc)
        return self.count_in_region(mp.map1, mp.map2)

    @staticmethod
    def get_filename():
        return "gap_in_kb_size.txt"



class KilobaseExperiment(Experiment):
    def __init__(self, map_size, sigma, pc, gap):
        Experiment.__init__(self, map_size, sigma, pc, gap)
        # self.gap *= math.sqrt(self.map_size) * self.sigma

    def out_border(self, map1, map2):
        alignment_ids = sorted(list(set(map1.ids) & set(map2.ids)))

        cut1 = map1.cuts
        cut2 = map2.cuts

        for i, id in enumerate(alignment_ids):
            if cut2[map2.id_nums[id]] - cut1[map1.id_nums[id]] - self.gap > 0:
                return True

            if cut2[map2.id_nums[id]] - cut1[map1.id_nums[id]] + self.gap < 0:
                return True

        return False


    def run(self):
        mp = MapManager(self.map_size, self.sigma, self.pc)
        return self.out_border(mp.map1, mp.map2)

    @staticmethod
    def get_filename():
        return "gap_in_kb.txt"


class Analyzer:

    def __init__(self, experiment, map_size, N):
        self.experimnet = experiment
        self.N = N
        self.map_size = map_size

    def analyze(self, gaps):
        pool = ThreadPool(4)
        results = pool.map(self.analyze_by_gap, gaps)
        pool.close()
        pool.join()
        self.save_results(results)

    def analyze_by_gap(self, gap):
        result = []
        for isigma in range(1, 10):
            sigma = isigma * 0.2
            for ipc in range(0, 8):
                pc = 0.6 + ipc * 0.05

                exp = self.experimnet(self.map_size, sigma, pc, gap)
                exp_result = 0
                for _ in range(self.N):
                    exp_result += exp.run()

                print("{} {} {} {}".format(sigma, pc, gap, exp_result / self.N))
                result.append("{} {} {} {}".format(sigma, pc, gap, exp_result / self.N))

        return result

    def optimal_gap(self, err_exp, size_exp, left, right, error_rate):
        result = []
        for isigma in range(1, 10):
            sigma = isigma * 0.2
            for ipc in range(0, 8):
                pc = 0.6 + ipc * 0.05

                exp = self.experimnet(err_exp, self.map_size, sigma, pc, error_rate)
                ## TODO SIGMA RIGHT
                opt_gap = exp.run(self.N, left, right)

                size_experimnet = size_exp(self.map_size, sigma, pc, opt_gap)
                cell_count = 0
                for _ in range(self.N):
                    cell_count += size_experimnet.run()

                print("{} {} {} {}".format(sigma, pc, round(opt_gap, 3), round(cell_count / self.N, 3)))
                result.append("{} {} {} {}".format(sigma, pc, round(opt_gap, 3), round(cell_count / self.N, 3)))
        self.save_result(result)

    def save_result(self, result):
        out = self.experimnet.get_filename()
        f = open(out, 'w')
        for r in result:
            print(r, file=f)
        f.close()

    def save_results(self, results):
        out = self.experimnet.get_filename()
        f = open(out, 'w')
        for output in results:
            for r in output:
                print(r, file=f)
        f.close()


class OptimalGapKB(Experiment):

    def __init__(self, experiment, map_size, sigma, pc, error_rate):
        Experiment.__init__(self, map_size, sigma, pc, 0)
        self.experiment = experiment
        self.error_rate = error_rate

    def binarySearch(self, left, right):
        gap = (left + right) / 2

        exp = self.experiment(self.map_size, self.sigma, self.pc, gap)
        errors = 0
        for _ in range(self.N):
            errors += exp.run()
        exp_result = errors / self.N

        if exp_result > self.error_rate:
             return self.binarySearch(gap, right)
        elif exp_result < self.error_rate:
            return self.binarySearch(left, gap)
        else:
            return gap

    def run(self, N, left, right):
        self.N = N
        return self.binarySearch(left, right)

    @staticmethod
    def get_filename():
        return "optimal_gap_in_kb.txt"


class OptimalGapFrag(Experiment):

    def __init__(self, experiment, map_size, sigma, pc, error_rate):
        Experiment.__init__(self, map_size, sigma, pc, 0)
        self.experiment = experiment
        self.error_rate = error_rate

    def binarySearch(self, left, right):

        for gap in range(left, self.map_size):
            exp = self.experiment(self.map_size, self.sigma, self.pc, gap)
            errors = 0
            for _ in range(self.N):
                errors += exp.run()
            exp_result = errors / self.N

            if exp_result <= self.error_rate:
                return gap


    def run(self, N, left, right):
        self.N = N
        return self.binarySearch(left, right)

    @staticmethod
    def get_filename():
        return "optimal_gap_in_frag.txt"

# exp = Experiment(25, 1.8, 0.6, 1, 0.5)
# map1 = Map([0, 1, 2, 5], [0, 1, 2, 5])
# map2 = Map([0, 1, 3, 4, 5], [0, 1, 3, 4, 5])
# print(exp.test_border_frags(map1, map2))
# print(exp.test_border_frags(map2, map1))

# exp = FragmentExperimentSize(1, 1, 1, 5)
# # print(exp.run_kb())
# map1 = Map([0, 1, 2, 5], [0, 1, 2, 5])
# map2 = Map([0, 1, 3, 4, 5], [0, 1, 3, 4, 5])
# print(exp.count_in_region(map1, map2))
# print(exp.count_in_region(map2, map1))

# ref_sum = 0
# res = []
# for r in reference_frags:
#     ref_sum += r
#     res.append(round(ref_sum, 3))
# print(res)
#

if __name__ == "__main__":

    # op = OptimalGap(KilobaseExperiment, 25, 0.6, 0.8, 0.01)
    # op.run(1000, 1, 20 * 0.8)
    # opa_kb = Analyzer(OptimalGapKB, 25, 100)
    # opa_kb.optimal_gap(KilobaseExperiment, KilobaseExperimentSize, 1, 40, 0.1)
    #
    # opa_frag = Analyzer(OptimalGapFrag, 25, 100)
    # opa_frag.optimal_gap(FragmentExperiment, FragmentExperimentSize, 1, 40, 0.1)

    fe = Analyzer(FragmentExperiment, 25, 100)
    fe.analyze([1, 2, 3, 4])
    # ke = Analyzer(KilobaseExperiment, 25, 100)
    # ke.analyze([1, 2, 3, 4])
    fes = Analyzer(FragmentExperimentSize, 25, 150)
    fes.analyze([1, 2, 3, 4])
    # kes = Analyzer(KilobaseExperimentSize, 25, 150)
    # kes.analyze([1, 2, 3, 4])


