from random import shuffle
from time import sleep

ranks = {'Ace': [1, 11], 'Two': 2, 'Three': 3, 'Four': 4, 'Five': 5, 'Six': 6, 'Seven': 7, 'Eight': 8, 'Nine': 9,
         'Ten': 10, 'Jack': 10, 'Queen': 10, 'King': 10}
suits = ['Clubs', 'Diamonds', 'Hearts', 'Spades']


class Card:
    def __init__(self, suit, rank):
        self.suit = suit
        self.rank = rank
        self.value = ranks[rank]

    def __str__(self):
        return self.rank + ' of ' + self.suit


class Deck:
    def __init__(self):
        used_deck = []
        complete_deck = []
        for suit in suits:
            for rank in ranks:
                new_card = Card(suit, rank)
                used_deck.append(new_card)
                complete_deck.append(new_card)

        self.current_deck = used_deck
        self.all_cards = complete_deck
        self.length = len(self.current_deck)

    def reset(self):
        self.current_deck = self.all_cards

    def mix(self):
        shuffle(self.all_cards)

    def deal(self):
        return self.all_cards.pop(0)


# Game-play Logic:
def black_jack():
    player_total = 100
    player_start = 100
    wins = 0
    losses = 0
    ties = 0

    # Introduction
    print('Hello, welcome to Python Casino, I will be your Black Jack dealer.\n')
    print(
        'Just a recap of the rules: Number cards are the value of their number, face cards are worth 10, and an ace is worth either 1 or 11.\nYou will see one dealer card then look at your hand and continue to hit or stay until satisfied with your total or you bust, if the dealer beats you or you bust, you lose your bet,\notherwise your bet is doubled.')

    game_on = True
    while game_on:
        # Set Up
        player_hand = []
        dealer_hand = []
        black_jack_deck = Deck()
        black_jack_deck.mix()

        # Formalities
        fair = ['']
        if player_total - player_start < 0:
            fair[0] = 'down'
        elif player_total - player_start >= 0:
            fair[0] = 'up'
        stat_check = True
        while stat_check:
            stat = input(
                '\nWould you like to know your statistics for today while I get the cards ready for the next game? ("yes" or "no"): ')
            if stat == 'yes':
                print(
                    f'You have won {wins} time(s) today, lost {losses} time(s), tied {ties} time(s), and are {fair[0]} ${abs(player_total - player_start)} today.')
                stat_check = False
                continue
            elif stat == 'no':
                stat_check = False
                continue
            else:
                print('Please enter either "Yes" or "No"')
                stat_check = True

        # Betting logistics
        print(f'\nI also see that you have ${player_total} in your account so place a bet you can afford.')
        bet_check = True
        while bet_check:
            char_bet_amount = input('Enter a betting amount: $')
            # Have to include the line that throws the error
            try:
                bet_amount = int(char_bet_amount)
                if int(char_bet_amount) > player_total:
                    print('Please lower the amount.')
                    bet_check = True
                elif int(char_bet_amount) <= 1:
                    print('Please enter a number greater than 0.')
                    bet_check = True
                elif player_total >= int(char_bet_amount) > 0:
                    print(
                        'Got it! Your bet will be added to the pool. If you win your earnings will be double your bet, if you lose the dealer will keep your bet.')
                    bet_check = False
            except ValueError:
                print('Please print an integer.')
                bet_check = True

        # Game Start
        # Dealer Set Up
        print("\nNow let's begin:")
        deal_one = black_jack_deck.deal()
        deal_two = black_jack_deck.deal()
        dealer_hand.append(deal_one)
        dealer_hand.append(deal_two)
        print("I'll go ahead and deal the cards.\n")
        sleep(4)

        # After being shown 1/2 dealer cards, user plays
        print(f"The dealers hand is an unknown card and the {dealer_hand[1]}.")
        # Deal user 2 cards, then ask if they want to hit or stay
        deal_three = black_jack_deck.deal()
        deal_four = black_jack_deck.deal()
        player_hand.append(deal_three)
        player_hand.append(deal_four)
        print(f"Your hand is {player_hand[0]}, and {player_hand[1]}.")

        # Ace Logistics out of initial hand
        if deal_three.rank == 'Ace':
            ace_check = True
            while ace_check:
                value_check = input(
                    "Since your first card is an ace, what do you want the value to be '1' or '11': ")
                if value_check == '1':
                    deal_four.value = 1
                    ace_check = False
                elif value_check == '11':
                    deal_four.value = 11
                    ace_check = False
                else:
                    print("Please enter the number '1' or '11'.")
                    ace_check = True
        elif deal_four.rank == 'Ace':
            ace_check = True
            while ace_check:
                value_check = input(
                    "Since your second card is an ace, what do you want the value to be '1' or '11': ")
                if value_check == '1':
                    deal_four.value = 1
                    ace_check = False
                elif value_check == '11':
                    deal_four.value = 11
                    ace_check = False
                else:
                    print("Please enter the number '1' or '11'.")
                    ace_check = True

        # Check total if stayed, and check if busted.
        hit = True
        while hit:
            action = input(f"\nBased on the cards you have, do you want to 'hit or 'stay': ")
            if action == 'stay':
                player_sum = 0
                for card in player_hand:
                    player_sum += card.value
                print(f'Your hand value: {player_sum}')
                hit = False
            elif action == 'hit':
                given_card = black_jack_deck.deal()
                player_hand.append(given_card)
                print(f'You drew a {given_card}!')
                # Logic to account for aces
                if given_card.rank == 'Ace':
                    ace_check = True
                    while ace_check:
                        value_check = input("Do you want your ace value to be '1' or '11': ")
                        if value_check == '1':
                            given_card.value = 1
                            ace_check = False
                        elif value_check == '11':
                            given_card.value = 11
                            ace_check = False
                        else:
                            print("Please enter the number '1' or '11'.")
                            ace_check = True
                # Back to logic for hitting
                player_sum = 0
                for card in player_hand:
                    player_sum += card.value
                print(f'Your hand value: {player_sum}')
                if player_sum > 21:
                    hit = False
                else:
                    hit = True
            else:
                print("Please use the right term ('hit or 'stay')")

        # If User didn't bust, then Dealer goes, showing the cards now, he goes until victory or bust.
        if player_sum > 21:
            print(f'Sorry you busted, meaning the dealer wins and you lose your bet amount of ${bet_amount}.')
            losses += 1
            player_total -= bet_amount
        elif player_sum <= 21:
            # Dealers turn
            print("\nNow I'll go: ")
            print(f"First I'll show you all my current cards: {dealer_hand[0]} and {dealer_hand[1]}")

            # Initial Ace Check
            dealer_sum = 0
            for card in dealer_hand:
                if card.rank == 'Ace':
                    card.value = 11
                    dealer_sum += 11
                else:
                    dealer_sum += card.value
            print(f'So my current deck total is {dealer_sum}.')

            # Conditionals to continue
            # Dealer deck value will never be less than player deck
            dealer_hit = True
            while dealer_hit:
                if dealer_sum < player_sum:
                    print('...just give me a second to draw.\n')
                    given_dealer_card = black_jack_deck.deal()
                    dealer_hand.append(given_dealer_card)
                    sleep(4)
                    print(f'Ok I drew a(n) {given_dealer_card}.')
                    if given_dealer_card.rank == 'Ace' and dealer_sum <= 10:
                        dealer_sum += 11
                        print(f'So now, {dealer_sum} is the value of my deck.')
                    elif given_dealer_card.rank == 'Ace' and dealer_sum > 10:
                        dealer_sum += 1
                        print(f'So now, {dealer_sum} is the value of my deck.')
                    else:
                        dealer_sum += given_dealer_card.value
                        print(f'So now, {dealer_sum} is the value of my deck.')
                elif dealer_sum >= player_sum:
                    dealer_hit = False

            # Check who won, then deal with bets.
            if dealer_sum > 21:
                print('\nI went over the limit of 21 so you take the win.')
                wins += 1
                player_total += bet_amount
            elif 21 >= dealer_sum > player_sum:
                print('\nCool, so I won that round.')
                print(f'Thanks for playing and since I won I will take your bet amount of ${bet_amount}')
                losses += 1
                player_total -= bet_amount
            elif player_sum > dealer_sum:
                print('\nCongratulations! You won the round!')
                winnings = bet_amount * 2
                print(f'Your winnings are ${winnings}')
                player_total += winnings
            elif player_sum == dealer_sum:
                print("\nOk I don't want to hit anymore, so the game will end in a tie meaning you keep your bet.")
                ties += 1

        # Close out the game by asking if they want to play again and reset game.
        if player_total > 0:
            play_game_checker = True
            while play_game_checker:
                cont = input('\nWould you like to play again? ("yes" or "no"): ')
                if cont == 'yes':
                    print('Ok awesome!')
                    black_jack_deck.reset()
                    # print('\n')
                    game_on = True
                    play_game_checker = False
                elif cont == 'no':
                    print('Ok, thanks for playing today!')
                    game_on = False
                    play_game_checker = False
                else:
                    print("Please enter either 'yes' or 'no'.")
                    play_game_checker = True
        else:
            print('\nYou have no money left to bet, please come again another day if you would like to play again.')
            game_on = False

black_jack()
