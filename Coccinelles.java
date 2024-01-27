import java.util.Random;
import java.util.Arrays;
import java.util.ArrayList;
import java.io.IOException;
import java.io.FileWriter;

/*
 * @author Richard Ho & Théo Phan
 */
class Coccinelles {

  public static void main(String[] args) {
    // Test Annexe
    System.out.println("Début du grand TestAnnexe");
    int[][] G = {
        { 2, 2, 3, 4, 2 },
        { 6, 5, 1, 2, 8 },
        { 1, 1, 10, 1, 1 },
    };
    afficherGrille(G, "G");
    System.out.print("Valeurs des chemins gloutons depuis les cases (0,d) : ");
    afficherNg(G);
    System.out.println("");

    System.out.println("Programmation dynamique, case de départ (0,0)");
    int[][][] MA = calculerMA(G, 0);
    int[][] M = MA[0];
    int[][] A = MA[1];
    afficherGrille(M, "M");
    System.out.print("Un chemin maximum : ");
    acnpm(M, A);
    System.out.println(" Valeur : " + optimal(G, 0));
    System.out.println("");

    System.out.println("Programmation dynamique, case de départ (0,1)");
    MA = calculerMA(G, 1);
    M = MA[0];
    A = MA[1];
    afficherGrille(M, "M");
    System.out.print("Un chemin maximum : ");
    acnpm(M, A);
    System.out.println(" Valeur : " + optimal(G, 1));
    System.out.println("");

    System.out.println("Chemins max. depuis toutes les cases de départ (0,d)");
    for (int d = 0; d < G[0].length; d++) {
      MA = calculerMA(G, d);
      M = MA[0];
      A = MA[1];
      System.out.print("Un chemin maximum : ");
      acnpm(M, A);
      System.out.println(" Valeur : " + optimal(G, d));
    }
    System.out.println("");

    int[] Ng = glouton(G);
    int[] Nmax = optimal(G);
    afficherTab(Ng, "Ng");
    afficherTab(Nmax, "Nmax");
    afficherGainsRelatifs(Nmax, Ng);
    System.out.println("");

    System.out.println("Début du grand Gains Relatifs /8");
    int[][] GRelatif = {
        { 2, 4, 3, 9, 6 },
        { 1, 10, 15, 1, 2 },
        { 2, 4, 11, 26, 66 },
        { 36, 34, 1, 13, 30 },
        { 46, 2, 8, 7, 15 },
        { 89, 27, 10, 12, 3 },
        { 1, 72, 3, 6, 6 },
        { 3, 1, 2, 4, 5 },
    };
    afficherGrille(GRelatif, "GRelatif");
    int[] Ng2 = glouton(GRelatif);
    int[] Nmax2 = optimal(GRelatif);
    afficherGainsRelatifs(Nmax2, Ng2);
    System.out.println("");

    // Validation Statistique
    System.out.println("VALIDATION STATISTIQUE");
    int n = 10000;
    System.out.println("nruns = " + n);
    int vMin = 5;
    int vMax = 16;
    System.out.println("L au hasard dans [" + vMin + ":" + vMax + "]");
    System.out.println("C au hasard dans [" + vMin + ":" + vMax + "]");
    System.out.println("Nb. de pucerons / case au hasard dans [0:L+C]");

    ArrayList<float[]> Gain = new ArrayList<float[]>();
    for (int i = 0; i < n; i++) {
      int[][] Grille = grilleAleatoire(vMin, vMax);
      int L = Grille.length - 1;
      int C = Grille[0].length - 1;
      float[] GainsRelatifs = gainRelatif(optimal(Grille), glouton(Grille));
      Gain.add(GainsRelatifs);
      if (i % 100 == 1) {
        System.out.printf("Run %d/%d, (L,C) = (%d,%d) \n", i, n, L, C);
      }

    }
    createCSV(Gain, "gainsRelatifs_nruns=10000_Linf=5_Lsup=16_Cinf=6_Csup=16");
    System.out.println("GAINS.lenght=" + gainsLength(Gain) + ", min=" + minimum(Gain) + ", max=" + maximum(Gain)
        + ", mean=" + moyenne(Gain) + ", med=" + mediane(Gain) + ", ecart-type=" + ecartType(Gain) + "\n");
    Gain.clear();

    ArrayList<float[]> Gain2 = new ArrayList<float[]>();
    System.out.println("");
    n = 100000;
    System.out.println("nruns = " + n);
    vMin = 2;
    vMax = 16;
    System.out.println("L au hasard dans [" + vMin + ":" + vMax + "]");
    System.out.println("C au hasard dans [" + vMin + ":" + vMax + "]");
    System.out.println("Nb. de pucerons / case au hasard dans [0:L*C]");

    for (int i = 0; i < n; i++) {
      int[][] Grille2 = grilleAleatoireComplexe(vMin, vMax);
      int L = Grille2.length - 1;
      int C = Grille2[0].length - 1;
      for (int l = 0; l < L + 1; l++) {
        permutationAleatoire(Grille2[l]);
      }
      float[] GainsRelatifs2 = gainRelatif(optimal(Grille2), glouton(Grille2));
      Gain2.add(GainsRelatifs2);
      if (i % 1000 == 1) {
        System.out.printf("Run %d/%d, (L,C) = (%d,%d) \n", i, n, L, C);
      }
    }
    createCSV(Gain2, "gainsRelatifs_nruns=100000_Linf=2_Lsup=16_Cinf=2_Csup=16");
    System.out.println("GAINS.lenght=" + gainsLength(Gain2) + ", min=" + minimum(Gain2) + ", max=" + maximum(Gain2)
        + ", mean=" + moyenne(Gain2) + ", med=" + mediane(Gain2) + ", ecart-type=" + ecartType(Gain2) + "\n");
    Gain2.clear();
  }

  // Outils

  static void afficherNg(int[][] tab) {
    int[] Ng = glouton(tab);
    System.out.println("Ng = " + Arrays.toString(Ng));
  }

  static void afficherNmax(int[][] tab) {
    int[] Nmax = optimal(tab);
    System.out.println("Nmax = " + Arrays.toString(Nmax));
  }

  static void afficherTab(int[] tab, String name) {
    System.out.println("" + name + " = " + Arrays.toString(tab));
  }

  static void afficherGainsRelatifs(int[] Nmax, int[] Ng) {
    float[] Gain = gainRelatif(Nmax, Ng);
    System.out.println("Gains relatifs = " + Arrays.toString(Gain));
  }

  static void afficherGrille(int[][] tab, String nom) {
    System.out.println("Grille " + nom + " :");
    for (int l = tab.length - 1; l >= 0; l--) {
      System.out.println("" + nom + "[" + l + "] : " + Arrays.toString(tab[l]));
    }
  }

  static int argMax(int[][] M) {
    int argMax = 0;
    int Lmax = M.length - 1;
    for (int c = 0; c < M[Lmax].length; c++) {
      if (M[Lmax][argMax] < M[Lmax][c]) {
        argMax = c;
      }
    }
    return argMax;
  }

  static void createCSV(ArrayList<float[]> Gain, String name) {
    try (FileWriter csvWriter = new FileWriter(name + ".csv")) {
      for (float[] GainsRelatifs : Gain) {
        String ligne = "";
        for (float GainRelatif : GainsRelatifs) {
          ligne += "" + GainRelatif + "\n";
        }
        csvWriter.append(ligne);
      }
      csvWriter.close();
      System.out.println("Les gains relatifs sont dans le fichier " + name + ".csv");
    } catch (IOException e) {
      System.out.println("Erreur de création du fichier CSV");
    }
  }

  // Naif Version
  static int glouton(int[][] G, int d) {
    int L = G.length;
    int C = G[0].length;
    int c = d;
    int puceM = G[0][c];

    for (int l = 0; l < L - 1; l++) {
      if (c == 0) {
        if (G[l + 1][c] < G[l + 1][c + 1]) {
          c++;
        }
      } else if (c == (C - 1)) {
        if (G[l + 1][c] < G[l + 1][c - 1]) {
          c--;
        }
      } else {
        if (G[l + 1][c] < G[l + 1][c - 1]) {
          if (G[l + 1][c - 1] < G[l + 1][c + 1]) {
            c++;
          } else {
            c--;
          }
        } else {
          if (G[l + 1][c] < G[l + 1][c + 1]) {
            c++;
          }
        }
      }
      puceM += G[l + 1][c];
    }
    return puceM;
  }

  static int[] glouton(int[][] G) {
    int C = G[0].length;
    int[] Ng = new int[C];
    for (int d = 0; d < C; d++) {
      Ng[d] = glouton(G, d);
    }
    return Ng;
  }

  /*
   * Dynamique
   *
   * Equation de récurrence des valeurs m(l,c) le nombre maximum de pucerons que
   * la coccinelle pourra manger à la ligne l et colonne c en ayant atterri sur
   * la case (0, d)
   * avec 0 <= l < L, 0 <= c < C et d une valeur fixe
   *
   * Soit v(l,c) la valeur/nombre de puceron à la ligne l et colonne c
   * avec 0 <= l < L et 0 <= c < C
   * 
   * - Base
   * m(0,c) = -1 , pour 0 <= c < C avec c != d
   * m(0,d) = v(0,d)
   * 
   * - Hérédité
   * m(l,c) = max( m(l-1,c') ) + v(l,c)
   * pour 1 <= l < L ; 0 <= c < C ; 0 <= c-1 <= c' <= c+1 < C
   * 
   */

  // C>1
  static int[][][] calculerMA(int[][] G, int d) {
    int L = G.length;
    int C = G[0].length;
    int[][] M = new int[L][C];
    int[][] A = new int[L][C];

    // Base
    for (int c = 0; c < C; c++) {
      M[0][c] = -1;
      A[0][c] = -1;
    }
    M[0][d] = G[0][d];
    A[0][d] = 0;

    // Heredite
    for (int l = 1; l < L; l++) {
      for (int c = 0; c < C; c++) {
        A[l][c] = -1;
        int alc;
        if (c == 0) {
          if (M[l - 1][c] > M[l - 1][c + 1]) {
            M[l][c] = M[l - 1][c];
            alc = c;
          } else {
            M[l][c] = M[l - 1][c + 1];
            alc = c + 1;
          }
          if (M[l][c] != -1) {
            M[l][c] += G[l][c];
            A[l][c] = alc;
          }
        } else if (c == (C - 1)) {
          if (M[l - 1][c] > M[l - 1][c - 1]) {
            M[l][c] = M[l - 1][c];
            alc = c;
          } else {
            M[l][c] = M[l - 1][c - 1];
            alc = c - 1;
          }
          if (M[l][c] != -1) {
            M[l][c] += G[l][c];
            A[l][c] = alc;
          }
        } else {
          if (M[l - 1][c] > M[l - 1][c - 1] && M[l - 1][c] > M[l - 1][c + 1]) {
            M[l][c] = M[l - 1][c];
            alc = c;
          } else if (M[l - 1][c - 1] > M[l - 1][c + 1]) {
            M[l][c] = M[l - 1][c - 1];
            alc = c - 1;
          } else {
            M[l][c] = M[l - 1][c + 1];
            alc = c + 1;
          }
          if (M[l][c] != -1) {
            M[l][c] += G[l][c];
            A[l][c] = alc;
          }
        }
      }
    }
    return new int[][][] { M, A };
  }

  static void acnpm(int[][] A, int l, int c) {
    if (l == 0) {
      System.out.print("(0," + c + ")");
      return;
    }
    int i = A[l][c];
    acnpm(A, l - 1, i);
    System.out.print("--->" + "(" + l + "," + c + ")");
  }

  static void acnpm(int[][] M, int[][] A) {
    int L = M.length;
    int cStar = argMax(M);
    acnpm(A, L - 1, cStar);
  }

  static int optimal(int[][] G, int d) {
    int[][][] MA = calculerMA(G, d);
    int[][] M = MA[0];
    int L = M.length;
    return M[L - 1][argMax(M)];
  }

  static int[] optimal(int[][] G) {
    int[] Nmax = new int[G[0].length];
    for (int d = 0; d < G[0].length; d++) {
      Nmax[d] = optimal(G, d);
    }
    return Nmax;
  }

  static float[] gainRelatif(int[] Nmax, int[] Ng) {
    int C = Nmax.length;
    float[] Gain = new float[C];

    for (int d = 0; d < C; d++) {
      Gain[d] = (float) (Nmax[d] - Ng[d]) / Ng[d];
    }
    return Gain;
  }

  /// Outils Statistique
  static int[] permutationAleatoire(int[] T) {
    int n = T.length;
    // Calcule dans T une permutation aleatoire de T et retourne T
    Random rand = new Random(); // bibliotheque java.util.Random
    for (int i = n; i > 0; i--) {
      int r = rand.nextInt(i); // r est au hasard dans [0:i]
      permuter(T, r, i - 1);
    }
    return T;
  }

  static void permuter(int[] T, int i, int j) {
    int ti = T[i];
    T[i] = T[j];
    T[j] = ti;
  }

  static int[][] grilleAleatoire() {
    Random random = new Random();
    int L = random.nextInt(11) + 5;
    int C = random.nextInt(11) + 5;
    int[][] G = new int[L][C];
    for (int l = 0; l < L; l++) {
      for (int c = 0; c < C; c++) {
        G[l][c] = random.nextInt(L + C - 1) + 1;
      }
    }
    return G;
  }

  static int[][] grilleAleatoire(int min, int max) {
    Random random = new Random();
    int L = random.nextInt(max - min) + min;
    int C = random.nextInt(max - min) + min;
    int[][] G = new int[L][C];
    for (int l = 0; l < L; l++) {
      for (int c = 0; c < C; c++) {
        G[l][c] = random.nextInt(L + C - 1) + 1;
      }
    }
    return G;
  }

  static int[][] grilleAleatoireComplexe(int min, int max) {
    Random random = new Random();
    int L = random.nextInt(max - min) + min;
    int C = random.nextInt(max - min) + min;
    int[][] G = new int[L][C];
    for (int l = 0; l < L; l++) {
      for (int c = 0; c < C; c++) {
        G[l][c] = random.nextInt(L * C - 1) + 1;
      }
    }
    return G;
  }

  static float minimum(ArrayList<float[]> Gain) {
    float min = Float.MAX_VALUE;
    for (int i = 0; i < Gain.size(); i++) {
      float[] gainRelatif = Gain.get(i);
      for (int j = 0; j < gainRelatif.length; j++) {
        if (gainRelatif[j] < min) {
          min = gainRelatif[j];
        }
      }
    }

    return min;
  }

  static float maximum(ArrayList<float[]> Gain) {
    float max = Float.MIN_VALUE;
    for (int i = 0; i < Gain.size(); i++) {
      float[] gainRelatif = Gain.get(i);
      for (int j = 0; j < gainRelatif.length; j++) {
        if (gainRelatif[j] > max) {
          max = gainRelatif[j];
        }
      }
    }

    return max;
  }

  static float moyenne(ArrayList<float[]> Gain) {
    float somme = (float) 0;
    float nbGainRelatif = (float) 0;
    for (int i = 0; i < Gain.size(); i++) {
      float[] gainRelatif = Gain.get(i);
      for (int j = 0; j < gainRelatif.length; j++) {
        somme += gainRelatif[j];
      }
      nbGainRelatif += gainRelatif.length;
    }
    return somme / nbGainRelatif;
  }

  static float mediane(ArrayList<float[]> Gain) {
    int tailleGRTotal = 0;
    for (int i = 0; i < Gain.size(); i++) {
      float[] gainRelatif = Gain.get(i);
      for (int j = 0; j < gainRelatif.length; j++) {
        tailleGRTotal++;
      }
    }
    float[] GRTotal = new float[tailleGRTotal];
    int indice = 0;
    for (int i = 0; i < Gain.size(); i++) {
      float[] gainRelatif = Gain.get(i);
      for (int j = 0; j < gainRelatif.length; j++) {
        GRTotal[indice] = gainRelatif[j];
        indice++;
      }
    }
    qs(GRTotal, 0, GRTotal.length);
    if (GRTotal.length % 2 == 0) {
      return GRTotal[GRTotal.length / 2];
    } else {
      return GRTotal[GRTotal.length / 2 - 1];
    }
  }

  static void qs(float[] GainTrie, int i, int j) {
    if (j - i <= 1) {
      return;
    }
    int k = partition(GainTrie, i, j);
    qs(GainTrie, i, k);
    qs(GainTrie, k + 1, j);
  }

  static int partition(float[] GainTrie, int i, int j) {
    int k = i;
    int jp = j;
    while (jp != k + 1) {
      if (GainTrie[k + 1] <= GainTrie[k]) {
        permuter(GainTrie, k + 1, k);
        k++;
      } else {
        permuter(GainTrie, k + 1, j - 1);
        jp--;
      }
    }
    return k;

  }

  static void permuter(float[] GainTrie, int i, int j) {
    float ti = GainTrie[i];
    GainTrie[i] = GainTrie[j];
    GainTrie[j] = ti;
  }

  static float ecartType(ArrayList<float[]> Gain) {
    float moy = moyenne(Gain);
    float somme = (float) 0;
    float nbGainRelatif = (float) 0;
    for (int i = 0; i < Gain.size(); i++) {
      float[] gainRelatif = Gain.get(i);
      for (int j = 0; j < gainRelatif.length; j++) {
        somme += Math.pow(gainRelatif[j] - moy, 2);
      }
      nbGainRelatif += gainRelatif.length;
    }
    return (float) Math.sqrt(somme / nbGainRelatif);
  }

  static int gainsLength(ArrayList<float[]> Gain) {
    int tailleGRTotal = 0;
    for (int i = 0; i < Gain.size(); i++) {
      tailleGRTotal += Gain.get(i).length;
    }
    return tailleGRTotal;

  }
}
