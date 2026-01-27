/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "verifier.hpp"

#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stack>
#include <string>

using namespace std;

namespace {

void HoodOf(pg::Game& game, int src, std::ostream& logger) {
  logger << "=> ";
  for (const int* curedge = game.outs(src); *curedge != -1;
       curedge++) {
    int to = *curedge;
    logger << game.label_vertex(to) << ",";
  }
  logger << std::endl;
}

bool
processSCC(pg::Game& game,
           std::vector<int>& res,
           std::ostream& logger,
           int* done,
           int v,
           int prio)
{
  int max_prio = -1;
  int scc_size = 0;

  for (auto it = res.rbegin(); it != res.rend(); it++) {
    int n = *it;
    if (game.priority(n) > max_prio)
      max_prio = game.priority(n);
    scc_size++;
    done[n] = prio; // mark as done at prio
    if (n == v)
      break;
  }

  bool cycles = scc_size > 1 or game.getStrategy(v) == v or
                (game.getStrategy(v) == -1 and game.has_edge(v, v));

  if (cycles && (max_prio & 1) == (prio & 1)) {
    /**
     * Found! Report.
     */
    logger << "max_prio " << max_prio << " and prio " << prio << std::endl;
    logger << "scc_size " << scc_size << std::endl;
    logger << "plays a loop? " << (game.getStrategy(v) == v) << std::endl;
    logger << "implicit loop? "
           << (game.getStrategy(v) == -1 and game.has_edge(v, v)) << std::endl;
    logger << "\033[1;31mscc where loser wins\033[m with priority \033[1;34m"
           << max_prio << "\033[m\n";
    for (auto it = res.rbegin(); it != res.rend(); it++) {
      int n = *it;
      logger << " " << game.label_vertex(n);
      HoodOf(game, n, logger);
      if (n == v)
        break;
    }
    logger << std::endl;
    return true;
  }
  return false;
}

void
tarjanSearch(pg::Game& game, bool even, bool odd, std::ostream& logger)
{
  const int n_vertices = game.vertexcount();

  // Allocate datastructures for Tarjan search
  int* done = new int[n_vertices];
  int64_t* low = new int64_t[n_vertices];
  for (int i = 0; i < n_vertices; i++) {
    done[i] = -1;
    low[i] = 0;
  }

  std::vector<int> res;
  std::stack<int> st;

  int64_t pre = 0;

  for (int top = n_vertices - 1; top >= 0; top--) {
    // only if a dominion
    if (!game.isSolved(top))
      continue;

    int prio = game.priority(top);
    int winner = game.getWinner(top);

    // only compute SCC for a (probably) top vertex
    if (winner == 0 and !even)
      continue; // don't check even dominions
    if (winner == 1 and !odd)
      continue; // don't check odd dominions

    // only try to find an SCC where the loser wins
    if (winner == (prio & 1))
      continue;

    // only run the check if not yet done at priority <prio>
    if (done[top] == prio)
      continue;

    // set <bot> (in tarjan search) to current pre
    int64_t bot = pre;

    // start the tarjan search at vertex <top>
    st.push(top);

    while (!st.empty()) {
      int v = st.top();

      /**
       * When we see it for the first item, we assign the next number to it
       * and add it to <res>.
       */
      if (low[v] <= bot) {
        if (std::numeric_limits<decltype(pre)>::max() == pre)
          throw std::runtime_error("Ran out of pre indices in Tarjan!");
        low[v] = ++pre;
        res.push_back(v);
      }

      /**
       * Now we check all outgoing (allowed) edges.
       * If seen earlier, then update "min"
       * If new, then 'recurse'
       */
      int min = low[v];
      bool pushed = false;
      int wrapped[2] = { game.getStrategy(v), -1 };
      for (const int* curedge = (game.getStrategy(v) != -1) ? wrapped
                                                            : game.outs(v);
           *curedge != -1;
           curedge++) {
        int to = *curedge;
        // skip if to higher priority
        if (game.priority(to) > prio)
          continue;
        // skip if already found scc (done[to] set to prio)
        if (done[to] == prio)
          continue;
        // check if visited in this search
        if (low[to] <= bot) {
          // not visited, add to <st> and break!
          // logger << " edge " << game.label_vertex(v) << " -> " <<
          // game.label_vertex(to) << std::endl;
          st.push(to);
          pushed = true;
          break;
        } else {
          // visited, update min
          if (low[to] < min)
            min = low[to];
        }
      }
      if (pushed)
        continue; // we pushed a new vertex to <st>...

      // there was no edge to a new vertex
      // check if we are the root of a SCC
      if (min < low[v]) {
        // not the root (outgoing edge to lower ranked vertex)
        low[v] = min;
        st.pop();
        continue;
      }

      /**
       * We're the root!
       * Now we need to figure out if we have cycles with p...
       * Also, mark every vertex in the SCC as "done @ search p"
       */
      if (processSCC(game, res, logger, done, v, prio)) {
        delete[] done;
        delete[] low;

        throw std::runtime_error("loser can win");
      }

      /**
       * No problem found
       */
      int n;
      do {
        n = res.back();
        res.pop_back();
      } while (v != n);
      st.pop();
    }
  }

  delete[] done;
  delete[] low;
}
} // namespace

namespace pg {

void
Verifier::verify(bool fullgame, bool even, bool odd)
{
  // ensure the vertices are ordered properly
  game.ensure_sorted();
  // ensure that the arrays are built
  game.build_in_array(false);

  const int n_vertices = game.vertexcount();

  /**
   * The first loop removes all edges from "won" vertices that are not the
   * strategy. This turns each dominion into a single player game. Also some
   * trivial checks are performed.
   */
  for (int v = 0; v < n_vertices; v++) {
    // (for full solutions) check whether every vertex is won
    if (!game.isSolved(v)) {
      if (fullgame)
        throw std::runtime_error("not every vertex is won");
      else
        continue;
    }

    const bool winner = game.getWinner(v);

    if (winner == 0 and !even)
      continue; // whatever
    if (winner == 1 and !odd)
      continue; // whatever

    if (winner == game.owner(v)) {
      // if winner, check whether the strategy stays in the dominion
      int str = game.getStrategy(v);
      if (str == -1) {
        throw std::runtime_error("winning vertex has no strategy");
      } else if (!game.has_edge(v, str)) {
        throw std::runtime_error("strategy is not a valid move");
      } else if (!game.isSolved(str) or game.getWinner(str) != winner) {
        throw std::runtime_error("strategy leaves dominion");
      }
      n_strategies++; // number of checked strategies
    } else {
      // if loser, check whether the loser can escape
      for (auto curedge = game.outs(v); *curedge != -1; curedge++) {
        int to = *curedge;
        if (!game.isSolved(to) or game.getWinner(to) != winner) {
          logger << "escape edge from " << game.label_vertex(v) << " to "
                 << game.label_vertex(to) << std::endl;
          throw std::runtime_error("loser can escape");
        }
      }
      // and of course check that no strategy is set
      if (game.getStrategy(v) != -1)
        throw std::runtime_error("losing vertex has strategy");
    }
  }

  tarjanSearch(game, even, odd, logger);
}

} // namespace pg
