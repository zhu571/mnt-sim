"""Phase 3 fragment recognition and de-excitation checks."""

from __future__ import annotations

import pathlib
import sys
import os

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mnt_sim.imqmd import (  # noqa: E402
    GaussianPacket,
    ImQMDNucleus,
    benchmark_evaporation_chain,
    coalescence_light,
    compute_fragment_excitation,
    grid_energy_diagnostics,
    evaporate_chain,
    fission_competition,
    impact_parameter_scan,
    initialize_nucleus,
    minimum_spanning_tree,
    propagate,
    reaction_fragments,
    validate_neutron_separation_energies,
    weisskopf_evaporation,
)


def _merge_nuclei(left: ImQMDNucleus, right: ImQMDNucleus, separation: float = 30.0) -> ImQMDNucleus:
    packets: list[GaussianPacket] = []
    for source, shift in ((left, -0.5 * separation), (right, 0.5 * separation)):
        for packet in source.packets:
            packets.append(
                GaussianPacket(
                    packet.r_i + np.array([shift, 0.0, 0.0]),
                    packet.p_i.copy(),
                    packet.sigma_r,
                    packet.is_proton,
                )
            )
    return ImQMDNucleus(left.Z + right.Z, left.N + right.N, packets, edf=left.edf)


def test_mst_two_fragments():
    left = initialize_nucleus(2, 4, sigma_r=1.1, seed=11)
    right = initialize_nucleus(2, 4, sigma_r=1.1, seed=12)
    system = _merge_nuclei(left, right, separation=25.0)

    fragments = minimum_spanning_tree(system, r_cut=3.0, p_cut=None)
    heavy = sorted(fragment.A for fragment in fragments if fragment.A > 1)
    assert heavy == [4, 4]


def test_coalescence_deuteron():
    packets = [
        GaussianPacket(np.array([0.0, 0.0, 0.0]), np.array([20.0, 0.0, 0.0]), 1.1, True),
        GaussianPacket(np.array([1.2, 0.0, 0.0]), np.array([30.0, 0.0, 0.0]), 1.1, False),
    ]
    nucleus = ImQMDNucleus(1, 1, packets)

    clusters = coalescence_light(nucleus)
    assert any(cluster.Z == 1 and cluster.A == 2 for cluster in clusters)


def test_cold_nucleus_excitation_is_zero():
    Z, A = 8, 16
    nucleus = initialize_nucleus(Z, A, sigma_r=1.1, seed=21)
    fragment = minimum_spanning_tree(nucleus, r_cut=3.0, p_cut=None)[0]

    assert abs(compute_fragment_excitation(nucleus, fragment)) < 0.1


def test_heated_nucleus_excitation_is_positive():
    Z, A = 8, 16
    nucleus = initialize_nucleus(Z, A, sigma_r=1.1, seed=21)
    fragment = minimum_spanning_tree(nucleus, r_cut=3.0, p_cut=None)[0]
    momenta = nucleus.momenta
    momenta[0, 0] += 50.0
    momenta[1, 0] -= 50.0
    nucleus.momenta = momenta

    assert compute_fragment_excitation(nucleus, fragment) > 0.0


def test_grid_energy_diagnostics_reports_scale():
    diagnostics = grid_energy_diagnostics(20, 40, sigma_r=1.1, seed=40)

    assert diagnostics["raw_total"] < diagnostics["target_total"]
    assert 0.35 <= diagnostics["nuclear_scale"] <= 1.25
    assert abs(diagnostics["scaled_total"] + diagnostics["energy_offset"] - diagnostics["target_total"]) < 5.0


def test_evaporation_reduces_A():
    channels = weisskopf_evaporation(20, 40, 20.0, n_max=6)
    evaporated = [(Z, A, p) for Z, A, p in channels if A < 40 and p > 0.0]
    mean_a = sum(A * p for _, A, p in channels)

    assert evaporated
    assert mean_a < 40.0


def test_evaporation_chain_cools_below_1mev():
    chain = evaporate_chain(28, 64, 20.0, rng=np.random.default_rng(7))

    assert chain.residue_Z <= 28
    assert chain.residue_A <= 64
    assert chain.residue_E_star < 1.0
    assert chain.steps


def test_fission_competition_rises_with_excitation():
    low = fission_competition(92, 238, 10.0)
    high = fission_competition(92, 238, 80.0)

    assert 0.0 <= low <= 1.0
    assert 0.0 <= high <= 1.0
    assert high >= low


def test_neutron_separation_matches_mass_table_sample():
    discrepancies = validate_neutron_separation_energies(
        nuclei=[(20, 40), (28, 64), (50, 132), (82, 208), (92, 238)],
        warn=False,
    )

    assert discrepancies == []


def test_evaporation_benchmark_reasonable():
    benchmark = benchmark_evaporation_chain(
        test_points=[(28, 64, 20.0), (92, 238, 60.0)],
        n_trials=32,
        seed=9,
    )
    by_nucleus = {(item["Z"], item["A"]): item for item in benchmark}

    ni = by_nucleus[(28, 64)]
    u = by_nucleus[(92, 238)]
    assert 0.0 <= ni["fission_fraction"] == 0.0
    assert ni["mean_residue_A"] <= 64.0
    assert ni["mean_steps"] >= 1.0
    assert 0.0 <= u["fission_fraction"] <= 1.0
    assert u["mean_residue_A"] <= 238.0


def test_cold_nucleus_few_fragments():
    nucleus = initialize_nucleus(20, 40, sigma_r=1.1, seed=40)
    fragments = minimum_spanning_tree(nucleus, r_cut=3.0, p_cut=None)

    assert len(fragments) == 1
    assert fragments[0].Z == 20
    assert fragments[0].A == 40


def test_two_nucleus_collision_separates():
    projectile = initialize_nucleus(8, 16, sigma_r=1.1, seed=161)
    target = initialize_nucleus(8, 16, sigma_r=1.1, seed=162)
    packets: list[GaussianPacket] = []
    separation = 24.0
    impact_parameter = 8.0
    beam_momentum = 220.0
    reference_positions: list[np.ndarray] = []
    reference_group_ids: list[int] = []

    for packet in projectile.packets:
        shift = np.array([-0.5 * separation, 0.5 * impact_parameter, 0.0])
        packets.append(
            GaussianPacket(
                packet.r_i + shift,
                packet.p_i + np.array([beam_momentum, 0.0, 0.0]),
                packet.sigma_r,
                packet.is_proton,
            )
        )
        reference_positions.append(packet.r_i + shift)
        reference_group_ids.append(0)
    for packet in target.packets:
        shift = np.array([0.5 * separation, -0.5 * impact_parameter, 0.0])
        packets.append(
            GaussianPacket(
                packet.r_i + shift,
                packet.p_i + np.array([-beam_momentum, 0.0, 0.0]),
                packet.sigma_r,
                packet.is_proton,
            )
        )
        reference_positions.append(packet.r_i + shift)
        reference_group_ids.append(1)

    system = ImQMDNucleus(16, 16, packets, edf=projectile.edf, reference_positions=np.asarray(reference_positions))
    system.reference_group_ids = np.asarray(reference_group_ids, dtype=int)
    propagate(
        system,
        dt=1.0,
        n_steps=150,
        sample_every=50,
        remove_cm_drift=False,
        use_surface_term=True,
        use_static_stabilizer=False,
    )
    fragments = reaction_fragments(system)
    heavy = [fragment for fragment in fragments if fragment.A >= 8]

    assert len(heavy) >= 2
    assert heavy[0].position[0] * heavy[1].position[0] < 0.0


def test_mnt_event_has_transfer():
    if os.environ.get("RUN_SLOW_MNT_TESTS") != "1":
        return

    projectile = initialize_nucleus(92, 238, sigma_r=1.1, seed=23801)
    target = initialize_nucleus(92, 238, sigma_r=1.1, seed=23802)
    packets: list[GaussianPacket] = []
    separation = 28.0
    impact_parameter = 4.0
    beam_momentum = np.sqrt(2.0 * 938.9 * 7.0)
    reference_positions: list[np.ndarray] = []
    reference_group_ids: list[int] = []

    for packet in projectile.packets:
        shift = np.array([-0.5 * separation, 0.5 * impact_parameter, 0.0])
        packets.append(
            GaussianPacket(
                packet.r_i + shift,
                packet.p_i + np.array([beam_momentum, 0.0, 0.0]),
                packet.sigma_r,
                packet.is_proton,
            )
        )
        reference_positions.append(packet.r_i + shift)
        reference_group_ids.append(0)
    for packet in target.packets:
        shift = np.array([0.5 * separation, -0.5 * impact_parameter, 0.0])
        packets.append(GaussianPacket(packet.r_i + shift, packet.p_i.copy(), packet.sigma_r, packet.is_proton))
        reference_positions.append(packet.r_i + shift)
        reference_group_ids.append(1)

    system = ImQMDNucleus(184, 292, packets, edf=projectile.edf, reference_positions=np.asarray(reference_positions))
    system.reference_group_ids = np.asarray(reference_group_ids, dtype=int)
    system._collision_rng = np.random.default_rng(23804)
    propagate(
        system,
        dt=1.0,
        n_steps=600,
        sample_every=200,
        with_collisions=True,
        collision_dt=1.0,
        remove_cm_drift=False,
        use_surface_term=True,
        use_static_stabilizer=False,
    )
    fragments = reaction_fragments(system)
    heavy = [fragment for fragment in fragments if fragment.A > 100]

    assert getattr(system, "collision_stats", {}).get("accepted", 0) > 0
    assert any(fragment.Z != 92 for fragment in heavy)


def test_krni_scan_pipeline_smoke():
    if os.environ.get("RUN_SLOW_MNT_TESTS") != "1":
        return

    result = impact_parameter_scan(
        projectile_z=36,
        projectile_a=86,
        target_z=28,
        target_a=64,
        energy_per_a=25.0,
        b_values=[2.0, 4.0, 6.0],
        events_per_b=10,
        n_workers=1,
        time_fm_c=400.0,
        dt=1.0,
        collision_dt=1.0,
        fragment_method="iso-mst",
    )

    assert len(result.impact_parameter_results) == 3
    assert result.total_cross_section > 0.0
    assert any(z not in {28, 36} for z in result.dsigma_dz)


if __name__ == "__main__":
    for test in (
        test_mst_two_fragments,
        test_coalescence_deuteron,
        test_cold_nucleus_excitation_is_zero,
        test_heated_nucleus_excitation_is_positive,
        test_evaporation_reduces_A,
        test_evaporation_chain_cools_below_1mev,
        test_fission_competition_rises_with_excitation,
        test_neutron_separation_matches_mass_table_sample,
        test_evaporation_benchmark_reasonable,
        test_cold_nucleus_few_fragments,
        test_two_nucleus_collision_separates,
        test_mnt_event_has_transfer,
        test_krni_scan_pipeline_smoke,
    ):
        test()
    print("Fragment and de-excitation checks passed")
