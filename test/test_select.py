from easytrajh5.select import parse_number_list


def test_it():
    for num_s, num_list in [
        ("1 2 3", [1, 2, 3]),
        ("50-100", list(range(50, 101))),
        ("2-3,5-7", [2, 3, 5, 6, 7]),
        ("7-9,11,13-15 19,21", [7, 8, 9, 11, 13, 14, 15, 19, 21]),
    ]:
        parse_list = parse_number_list(num_s)
        parse_list.sort()
        num_list.sort()
        assert parse_list == num_list
